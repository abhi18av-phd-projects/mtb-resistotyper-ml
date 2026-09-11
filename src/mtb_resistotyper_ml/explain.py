"""Turn a prediction into per-mutation reasoning.

Nothing here is a new model. The explanation is the model: because the deployed
class is a scaler plus a coefficient vector, each fired mutation's signed logit
contribution is exactly what moved the prediction, and the contributions sum to
the logit. What this module adds is the labelling a coefficient cannot carry --
whether a gene is a known determinant of *this* drug, or a co-selected passenger
riding along on multi-drug linkage.

The confidence statement is read from the bundle's ``model_card.json`` and is
never hardcoded here. A hardcoded tier table goes stale the moment a model is
retrained, and would then misreport the operating range of an artefact that
correctly declares its own.
"""

from __future__ import annotations

# Canonical determinants, gene -> drugs. Used ONLY to label on- versus off-target.
# The prediction never consults it.
DETERMINANTS: dict[str, set[str]] = {
    "RIF": {"rpoB"},
    "INH": {"katG", "inhA", "fabG1", "ahpC"},
    "EMB": {"embB", "embA", "embC", "embR"},
    "PZA": {"pncA", "rpsA", "panD"},
    "MXF": {"gyrA", "gyrB"},
    "LEV": {"gyrA", "gyrB"},
    "AMI": {"rrs", "eis"},
    "KAN": {"rrs", "eis"},
    "STM": {"rpsL", "rrs", "gid"},
    "ETH": {"ethA", "ethR", "inhA", "fabG1"},
    "LZD": {"rrl", "rplC"},
    "BDQ": {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"},
    "CFZ": {"Rv0678", "pepQ", "mmpL5", "mmpS5"},
    "DLM": {"ddn", "fbiA", "fbiB", "fbiC", "fgd1", "Rv0678"},
}


def split_feature(feature: str) -> tuple[str, str]:
    """``raw__<gene>_<garc>`` -> ``(gene, garc)``. GARC may contain '_'."""
    gene, _, garc = feature[len("raw__"):].partition("_")
    return gene, garc


def _role(on_target: bool, contribution: float) -> str:
    if on_target:
        return "primary driver"
    if abs(contribution) >= 0.2:
        return "co-selected passenger (multi-drug linkage)"
    return "weak contributor"


def explain(scored: dict, drug: str) -> dict:
    """Annotate a scored result with per-mutation roles and the operating range."""
    determinants = DETERMINANTS.get(drug, set())
    reasons = []
    for c in scored["contributions"]:
        if not c["feature"].startswith("raw__") or c["value"] == 0:
            continue
        gene, garc = split_feature(c["feature"])
        on_target = gene in determinants
        reasons.append({
            "mutation": f"{gene}@{garc}",
            "gene": gene,
            "coefficient": c["coef"],
            "logit_contribution": c["logit_contribution"],
            "direction": "resistance" if c["coef"] > 0 else "protective",
            "on_target": on_target,
            "role": _role(on_target, c["logit_contribution"]),
        })

    rng = scored.get("operating_range", {})
    tier = rng.get("tier")
    caveat = (
        f"{drug} sits in the '{tier}' tier at AUC {rng['auc']:.3f} under "
        f"{rng.get('metric', 'the declared metric')}. This is a Layer 2 call: it "
        f"is intended for variants the curated catalogue grades Unknown or Fail, "
        f"and never overrides a catalogue call."
        if tier and rng.get("auc") is not None else
        "This model bundle declares no operating range; treat the call as unvalidated."
    )
    return {
        "primary_drivers": [r["mutation"] for r in reasons
                            if r["role"] == "primary driver"],
        "passengers_flagged": [r["mutation"] for r in reasons
                               if "passenger" in r["role"]],
        "confidence": {"tier": tier, "auc": rng.get("auc"), "caveat": caveat},
        "reasons": reasons,
    }


def causal_chain(result: dict, card: dict | None = None) -> list[dict]:
    """The chain of reasoning from an observed variant to a model call.

    A probability and a list of coefficients is not a reason; it is the arithmetic
    a reason would have to survive. This assembles the steps that actually
    connect the two, and it is careful about which links are causal claims and
    which are only associations.

    The honest distinction is at step 4. An on-target contributor is a mutation
    in a gene with an established mechanism for this drug, and the model's weight
    on it is consistent with that mechanism. An off-target contributor carries
    weight because it travels with resistance in the training population, through
    co-selection and lineage structure, not because it causes it. Both move the
    prediction; only the first supports a mechanistic reading, and a chain that
    presented them alike would be laundering correlation into explanation.

    Step 3 states what the feature's presence in the schema does and does not
    establish. Each bundle records its selection procedure at drug level in
    ``model_card.json``; the per-feature votes are not in the artefact, so the
    chain says the mutation entered through that procedure rather than claiming
    a vote tally it cannot see.
    """
    card = card or {}
    training = card.get("training", {})
    rng = card.get("operating_range", {})
    prov = card.get("provenance", {})
    data = card.get("data", {})
    drug = result.get("drug")
    reasons = result.get("reasons", [])
    on = [r for r in reasons if r["on_target"]]
    off = [r for r in reasons if not r["on_target"]]

    chain: list[dict] = []

    chain.append({
        "step": 1, "stage": "observation",
        "claim": f"{len(reasons)} modelled mutation(s) were called in this isolate"
                 if reasons else "no modelled mutation was called in this isolate",
        "detail": [r["mutation"] for r in reasons],
        "basis": "variant calls named in GARC against NC_000962.3",
        "kind": "measurement",
    })

    cat = result.get("provenance", {}).get("catalogue") or {}
    chain.append({
        "step": 2, "stage": "layer 1 declined",
        "claim": f"the curated catalogue returned no grade for {drug}, so the question "
                 f"passed to the model"
                 if result.get("provenance", {}).get("catalogue_consulted")
                 else "the catalogue was not consulted, so this call rests on the model alone",
        "detail": {"catalogue": cat.get("catalogue"), "version": cat.get("version")},
        "basis": "piezo returned U, F, or no matching rule for this drug",
        "kind": "procedure",
    })

    chain.append({
        "step": 3, "stage": "admission to the model",
        "claim": "these mutations are features of this drug's model, having entered "
                 "through its concordance selection",
        "detail": {
            "selection": training.get("feature_selection"),
            "candidates_considered": training.get("n_features"),
            "features_retained": training.get("n_concordant"),
            "per_feature_votes": "not recorded in the bundle; selection is documented "
                                 "at drug level, so presence here shows the mutation "
                                 "survived that procedure, not which votes it won",
        },
        "basis": f"model_card.json training block for {drug}",
        "kind": "procedure",
    })

    chain.append({
        "step": 4, "stage": "contribution",
        "claim": "each mutation moves the log-odds by its standardised value times "
                 "its coefficient; these sum to the logit exactly, so the "
                 "contributions are the model rather than an approximation of it",
        "detail": {
            "causal_reading_supported": [
                {"mutation": r["mutation"], "gene": r["gene"],
                 "logit": r["logit_contribution"], "direction": r["direction"],
                 "why": f"{r['gene']} is an established determinant for {drug}"}
                for r in on],
            "association_only": [
                {"mutation": r["mutation"], "gene": r["gene"],
                 "logit": r["logit_contribution"], "role": r["role"],
                 "why": f"{r['gene']} has no established mechanism for {drug}; the weight "
                        f"reflects co-occurrence with resistance in the training "
                        f"population, not causation"}
                for r in off],
        },
        "basis": "the fitted coefficient vector",
        "kind": "mixed",
    })

    chain.append({
        "step": 5, "stage": "aggregation",
        "claim": f"contributions and the intercept give P(resistant) = "
                 f"{result.get('p_resistant')}",
        "detail": {"threshold": result.get("threshold"),
                   "call": result.get("prediction")},
        "basis": "sigmoid of the summed logit",
        "kind": "arithmetic",
    })

    band = result.get("reliability", {})
    chain.append({
        "step": 6, "stage": "bounding",
        "claim": f"the call is reported at reliability {band.get('band')}/"
                 f"{band.get('of')} within the model's measured operating range",
        "detail": {
            "auc": rng.get("auc"), "metric": rng.get("primary_metric"),
            "tier": rng.get("tier"),
            "capped_by_tier": band.get("capped_by_tier"),
            "capped_by_evidence": band.get("capped_by_evidence"),
            "fitted_on": {"source": data.get("source"), "release": data.get("release"),
                          "n_isolates": data.get("n_isolates")},
            "catalogue_edition_at_training": prov.get("catalogue_edition"),
        },
        "basis": "model_card.json operating_range",
        "kind": "qualifier",
    })

    if not on and off:
        chain.append({
            "step": 7, "stage": "warning",
            "claim": "no mutation with an established mechanism for this drug was "
                     "carried; the call rests entirely on associations",
            "detail": [r["mutation"] for r in off],
            "basis": "no contributor lies in a canonical determinant gene",
            "kind": "qualifier",
        })
    elif not reasons:
        chain.append({
            "step": 7, "stage": "warning",
            "claim": "no modelled mutation was carried, so the call is the model's "
                     "base rate and cannot distinguish an isolate with no resistance "
                     "mechanism from one carrying a mechanism this model has no "
                     "feature for",
            "detail": {"n_model_features": training.get("n_concordant")},
            "basis": "the prediction rests on the intercept alone",
            "kind": "qualifier",
        })
    return chain


def render_chain(chain: list[dict]) -> str:
    """The chain as text, for the report and the command line."""
    lines = []
    for s in chain:
        lines.append(f"  {s['step']}. {s['stage'].upper()}  [{s['kind']}]")
        lines.append(f"     {s['claim']}")
        d = s.get("detail")
        if s["stage"] == "contribution":
            for c in d["causal_reading_supported"]:
                lines.append(f"       + {c['mutation']:<22} logit {c['logit']:+.3f}   "
                             f"mechanism: {c['why']}")
            for c in d["association_only"]:
                lines.append(f"       ~ {c['mutation']:<22} logit {c['logit']:+.3f}   "
                             f"association only: {c['why']}")
        elif isinstance(d, list) and d:
            lines.append(f"       {', '.join(str(x) for x in d)}")
        elif isinstance(d, dict):
            for k, v in d.items():
                if v not in (None, "", [], {}):
                    lines.append(f"       {k}: {v}")
        lines.append(f"     basis: {s['basis']}")
    return "\n".join(lines)
