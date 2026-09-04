#!/usr/bin/env python3
"""Score the catalogue-unclassified cohort with Layer 2, and grade the answers.

Layer 1 is deliberately NOT re-run. The cohort is, by construction, the set the
catalogue returned U or F on; asking it again would return U or F again and cost
an hour of piezo. What is being measured here is only what Layer 2 adds where
Layer 1 stopped.

The output keeps every individual call rather than only the summary, because the
interesting rows are the ones a summary hides: an isolate the catalogue could not
grade, the model called resistant, and the MIC confirms it.
"""
import json, sys
from pathlib import Path

sys.path.insert(0, "/srv")
from mtb_resistotyper_ml import ModelBundle, resolve

COHORT = Path(sys.argv[1]); MODELS = Path(sys.argv[2]); OUT = Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)

summary, all_rows = {}, []
for f in sorted(COHORT.glob("*.json")):
    if f.name == "cohort_manifest.json":
        continue
    drug = f.stem
    bundle = ModelBundle.load(MODELS / drug)
    rows = json.loads(f.read_text())
    tp = fp = tn = fn = 0
    for r in rows:
        res = resolve(r["instance"], bundle, use_catalogue=False)
        call, truth = res["prediction"], r["phenotype"]
        tp += call == "R" and truth == "R"
        fp += call == "R" and truth == "S"
        tn += call == "S" and truth == "S"
        fn += call == "S" and truth == "R"
        all_rows.append({
            "sample_id": r["sample_id"], "drug": drug,
            "catalogue_call": r["catalogue_call"], "phenotype": truth,
            "mic": r["mic"], "lineage": r["lineage"],
            "model_call": call, "p_resistant": res["p_resistant"],
            "tier": res["confidence"]["tier"], "auc": res["confidence"]["auc"],
            "reliability": res["reliability"]["band"],
            "n_carried": res["n_mutations_carried"],
            "drivers": res["primary_drivers"],
            "passengers": res["passengers_flagged"],
            "reasons": [{"mutation": x["mutation"], "logit": x["logit_contribution"],
                         "on_target": x["on_target"], "role": x["role"]}
                        for x in res["reasons"]],
            "outcome": ("true_positive" if call == "R" and truth == "R" else
                        "false_positive" if call == "R" else
                        "false_negative" if truth == "R" else "true_negative"),
        })
    n_r, n_s = tp + fn, tn + fp
    summary[drug] = {
        "n": len(rows), "n_resistant": n_r, "n_susceptible": n_s,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        # Sensitivity is the number that matters here: these are isolates the
        # catalogue could not grade, so a missed resistant call is a resistant
        # isolate nothing in the stack flagged.
        "sensitivity": round(tp / n_r, 4) if n_r else None,
        "specificity": round(tn / n_s, 4) if n_s else None,
        "ppv": round(tp / (tp + fp), 4) if (tp + fp) else None,
        "tier": ModelBundle.load(MODELS / drug).operating_range.get("tier"),
    }
    s = summary[drug]
    print(f"  {drug:5} n={len(rows):5,} R={n_r:4,}  sens={s['sensitivity']}  "
          f"spec={s['specificity']}  ppv={s['ppv']}  (tp={tp} fn={fn} fp={fp})", flush=True)

(OUT / "predictions.json").write_text(json.dumps(all_rows))
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))
tot_tp = sum(v["tp"] for v in summary.values()); tot_r = sum(v["n_resistant"] for v in summary.values())
print(f"OVERALL recovered {tot_tp}/{tot_r} resistant isolates the catalogue could not grade "
      f"({100*tot_tp/tot_r:.1f}%)", flush=True)
