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

from typing import Any

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
