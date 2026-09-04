"""Where the WHO catalogue stops: what Layer 2 adds, and what it costs.

The two-layer design claims the models earn their place by answering for variants
the curated catalogue cannot grade. Until now that was an assertion. This
dashboard tests it on the only cohort that can: the CRyPTIC v3.4.0 isolates the
WHO catalogue returned Unknown or Fail on, restricted to those carrying a
measured MIC, so every model call has ground truth beside it.

The finding it exists to show is not the one that was expected. At the threshold
the models ship with, Layer 2 recovers 6 of 862 resistant isolates. The ranking
is not the problem -- AUC on this cohort reaches 0.83 for rifampicin -- the
operating point is. The models were fitted where 38.5% of isolates are resistant
and deployed where 11% are, so their base rate sits far below 0.5 and almost
nothing crosses it. Recalibrating on this cohort, holding specificity at 0.90,
recovers 137 instead of 6.

That is a defect in the artefact contract rather than in the models: the card
declares an operating RANGE and no operating POINT, and the operating point turns
out to be what decides whether the second layer does anything at all.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

DATA = Path(os.environ.get("MTB_EVAL_DATA", "/data"))

st.set_page_config(page_title="Layer 2 vs the catalogue's blind spot",
                   layout="wide", initial_sidebar_state="expanded")

# Status palette, used only for status, and never without a label beside it.
GOOD, CRIT, WARN, CAT, MOD = "#0ca30c", "#d03b3b", "#fab219", "#2a78d6", "#1baf7a"


@st.cache_data
def load():
    summary = json.loads((DATA / "summary.json").read_text())
    preds = pd.DataFrame(json.loads((DATA / "predictions.json").read_text()))
    cohort = json.loads((DATA / "cohort_manifest.json").read_text()) \
        if (DATA / "cohort_manifest.json").exists() else {}
    return summary, preds, cohort


summary, preds, cohort = load()
drugs = sorted(summary)

st.title("Where the WHO catalogue stops")
st.caption("CRyPTIC v3.4.0 · WHO-UCN-GTB-PCI-2023.5 · isolates graded Unknown or Fail, "
           "restricted to those with a measured MIC")

st.warning("**Research use only. Not a diagnostic.** These models have not been "
           "prospectively validated and have no regulatory clearance. Nothing on this "
           "page should guide the treatment of a patient.", icon="⚠️")

tot_n = sum(v["n"] for v in summary.values())
tot_r = sum(v["n_resistant"] for v in summary.values())
tot_tp = sum(v["tp"] for v in summary.values())
recal_tp = sum((v.get("recalibrated") or {}).get("tp", 0) for v in summary.values())

c1, c2, c3, c4 = st.columns(4)
c1.metric("Unclassified isolates with a phenotype", f"{tot_n:,}")
c2.metric("Phenotypically resistant", f"{tot_r:,}",
          help="The catalogue could not grade these and the organism is resistant. "
               "This is the population Layer 2 exists for.")
c3.metric("Recovered at the shipped threshold", f"{tot_tp}",
          delta=f"{100*tot_tp/tot_r:.1f}% of them", delta_color="off")
c4.metric("Recovered after recalibration", f"{recal_tp}",
          delta=f"+{recal_tp - tot_tp} at specificity 0.90")

tab1, tab2, tab3, tab4 = st.tabs(
    ["The gap", "Does Layer 2 close it?", "The threshold", "Individual calls"])

with tab1:
    st.subheader("How much the catalogue leaves ungraded")
    st.write(
        "Of 810,615 catalogue calls in CRyPTIC v3.4.0, 58,269 are Unknown and 9,875 Fail — "
        "8.4% ungraded. The table below is the part of that with a measured MIC, which is "
        "the only part on which a model can be checked rather than merely run.")
    gap = pd.DataFrame([{
        "drug": d, "isolates": v["n"], "resistant": v["n_resistant"],
        "susceptible": v["n_susceptible"],
        "prevalence": v.get("prevalence"),
        "carry a modelled mutation": cohort.get(d, {}).get("n_with_modelled_mutation"),
        "model features": cohort.get(d, {}).get("n_features"),
    } for d, v in summary.items()]).sort_values("resistant", ascending=False)
    st.dataframe(gap, hide_index=True, use_container_width=True,
                 column_config={"prevalence": st.column_config.NumberColumn(format="%.3f")})
    st.caption("Prevalence here is 11% overall. The models were fitted at 38.5%. That "
               "difference is the whole of tab three.")

with tab2:
    st.subheader("Layer 2 against the phenotype, at the threshold it ships with")
    perf = pd.DataFrame([{
        "drug": d, "n": v["n"], "resistant": v["n_resistant"],
        "sensitivity": v["sensitivity"], "specificity": v["specificity"],
        "TP": v["tp"], "FN": v["fn"], "FP": v["fp"],
        "AUC here": v.get("auc_on_unclassified"),
        "AUC on training cohort": v.get("auc"),
        "tier": v.get("tier"),
    } for d, v in summary.items()]).sort_values("AUC here", ascending=False)
    st.dataframe(perf, hide_index=True, use_container_width=True,
                 column_config={c: st.column_config.NumberColumn(format="%.4f")
                                for c in ("sensitivity", "specificity", "AUC here",
                                          "AUC on training cohort")})
    st.error(
        "**Sensitivity is near zero and specificity is near one: the models call almost "
        "everything susceptible on this cohort.** Read alone, that says the second layer "
        "adds nothing. The AUC column says otherwise — a ranking at 0.83 for RIF and 0.79 "
        "for EMB is not a model with nothing to say. The calls are wrong because the cut "
        "is in the wrong place, which is a different problem with a different fix.", icon="🔍")
    st.bar_chart(perf.set_index("drug")[["AUC here"]], height=260,
                 color=MOD, y_label="AUC on catalogue-unclassified isolates")
    st.caption("0.5 is the no-information line. AMI (0.55) and KAN (0.51) sit on it — for "
               "those two the models genuinely have nothing to add here, and no threshold "
               "rescues them.")

with tab3:
    st.subheader("The operating point, and what the shipped one costs")
    d = st.selectbox("Drug", drugs, index=drugs.index("RIF") if "RIF" in drugs else 0)
    v = summary[d]
    roc = pd.DataFrame(v["roc"])
    rec = v.get("recalibrated")

    a, b, c = st.columns(3)
    a.metric("AUC on this cohort", v.get("auc_on_unclassified"))
    b.metric("Sensitivity at p ≥ 0.5", f"{v['sensitivity']:.3f}" if v["sensitivity"] is not None else "—",
             delta=f"{v['tp']} of {v['n_resistant']} resistant", delta_color="off")
    if rec:
        c.metric(f"Sensitivity at p ≥ {rec['threshold']:.4f}", f"{rec['sensitivity']:.3f}",
                 delta=f"{rec['tp']} of {v['n_resistant']} · specificity {rec['specificity']:.3f}")

    st.line_chart(roc.set_index("threshold")[["sensitivity", "specificity"]],
                  height=300, color=[CRIT, CAT])
    st.caption(
        "Sensitivity and specificity against the decision threshold. The shipped cut is "
        "0.5, far to the right of anywhere useful on this cohort: the model's base rate "
        "here sits below it, so nearly every isolate is called susceptible regardless of "
        "what it carries.")
    st.info(
        "**Why this is a contract defect, not a tuning detail.** `model_card.json` declares "
        "an operating RANGE — an AUC and a tier — and no operating POINT. A consumer "
        "therefore inherits 0.5 by default, which is only right where prevalence matches "
        "the training cohort. The card should carry a threshold, the prevalence it assumes, "
        "and the sensitivity and specificity it buys.", icon="📋")

with tab4:
    st.subheader("The isolates the catalogue could not grade")
    f1, f2, f3 = st.columns(3)
    dsel = f1.multiselect("Drug", drugs, default=drugs)
    osel = f2.multiselect("Outcome", sorted(preds["outcome"].unique()),
                          default=["true_positive", "false_negative"])
    only_carrying = f3.checkbox("Only isolates carrying a modelled mutation", value=True)

    view = preds[preds["drug"].isin(dsel) & preds["outcome"].isin(osel)]
    if only_carrying:
        view = view[view["n_carried"] > 0]
    st.write(f"{len(view):,} isolates")
    st.dataframe(
        view[["sample_id", "drug", "catalogue_call", "phenotype", "mic", "model_call",
              "p_resistant", "reliability", "lineage", "n_carried", "drivers", "outcome"]]
        .sort_values("p_resistant", ascending=False),
        hide_index=True, use_container_width=True, height=380)

    st.markdown("##### Why a call was made")
    if len(view):
        pick = st.selectbox("Isolate", view["sample_id"].unique()[:400])
        for r in view[view["sample_id"] == pick].to_dict("records"):
            st.markdown(f"**{r['drug']}** — catalogue `{r['catalogue_call']}` · "
                        f"phenotype `{r['phenotype']}` (MIC {r['mic']}) · "
                        f"model `{r['model_call']}` at P(R) {r['p_resistant']}")
            if r["reasons"]:
                st.dataframe(pd.DataFrame(r["reasons"]), hide_index=True,
                             use_container_width=True)
            else:
                st.caption("No modelled mutation carried; the call rests on the intercept "
                           "alone, and cannot distinguish an isolate with no resistance "
                           "mechanism from one carrying a mechanism this model has no "
                           "feature for.")

st.divider()
st.caption(
    "Layer 1 is the curated WHO catalogue read by piezo; Layer 2 is a versioned model "
    "bundle. Layer 1 is not re-run here — this cohort is by construction the set it "
    "returned Unknown or Fail on. Reproduce any row with "
    "`pip install mtb-resistotyper-ml`.")
