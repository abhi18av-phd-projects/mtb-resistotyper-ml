"""Scoring for a deployed per-drug model bundle.

The deployed model is a StandardScaler followed by an L2-logistic regression, so
scoring is a closed form:

    p(R) = sigmoid( ((x - scaler_mean) / scaler_scale) . coef + intercept )

Two consequences follow, and both are the reason this model class was chosen.
The whole artefact is a few kilobytes, and the per-feature logit contributions
(coef_i * z_i) are an *exact* attribution rather than a post-hoc approximation of
one -- prediction and explanation come from the same expression.

This module deliberately imports nothing beyond the standard library, so that the
scoring path can be audited, vendored, or reimplemented in another language
without pulling in a numerical stack.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


class ModelBundle:
    """A model.json plus the model_card.json that states where it works."""

    def __init__(self, model: dict, card: dict | None = None,
                 schema: dict | None = None) -> None:
        self.model, self.card, self.schema = model, card or {}, schema or {}

    @classmethod
    def load(cls, directory: str | Path) -> ModelBundle:
        d = Path(directory)
        def _read(name: str) -> dict | None:
            p = d / name
            return json.loads(p.read_text()) if p.exists() else None
        model = _read("model.json")
        if model is None:
            raise FileNotFoundError(f"no model.json in {d}")
        return cls(model, _read("model_card.json"), _read("feature_schema.json"))

    @property
    def drug(self) -> str:
        return self.model["drug"]

    @property
    def operating_range(self) -> dict:
        """The measured range this model may be applied over.

        Read from the card rather than hardcoded, so a retrained model that ships
        a different range is reported with ITS range and not with a stale one.
        """
        return self.card.get("operating_range", {})


# Reliability bands, after the reliability index of Ghosh et al.'s standalone
# MTB predictors (mcdr-mtb-standalone, bdqr-mtb-standalone), which map the
# margin between the top two class probabilities onto 1 to 5. The idea is worth
# taking: a laboratory reads a coarse band more reliably than it reads a
# probability, and the band costs nothing to compute.
#
# It is taken with one change. A margin measures how far apart this model held
# the two classes on this isolate. It says nothing about whether the model is
# any good, so an indifferent model can report a wide margin and look certain.
# The band is therefore capped by the tier the bundle declares, and a bundle
# that declares no range cannot report better than 2. The margin is separation,
# not calibration, and the cap is what keeps the two from being confused.
TIER_CEILING: dict[str, int] = {"usable": 5, "moderate": 4, "weak": 3}
UNDECLARED_CEILING = 2
# A call on which no modelled mutation fired rests on the intercept alone. It
# cannot distinguish "this isolate carries no resistance mechanism" from "this
# isolate carries one for which the model has no feature", and these bundles
# carry as few as seven features. Absence of evidence over a small feature set
# is weaker than presence of it, so such a call cannot reach the top band.
NO_EVIDENCE_CEILING = 4


def reliability(p_resistant: float, tier: str | None = None,
                n_fired: int | None = None) -> dict:
    """A 1 to 5 confidence band for one call, capped by what supports it."""
    margin = abs(2.0 * p_resistant - 1.0)
    raw = max(1, min(5, 1 if margin < 0.1 else math.ceil(round(margin * 10) / 2)))
    tier_ceiling = TIER_CEILING.get(tier or "", UNDECLARED_CEILING)
    no_evidence = n_fired == 0
    band = min(raw, tier_ceiling, NO_EVIDENCE_CEILING if no_evidence else 5)
    return {
        "band": band,
        "of": 5,
        "margin": round(margin, 4),
        "capped_by_tier": band < raw and tier_ceiling <= band,
        "capped_by_evidence": no_evidence and band < raw,
        "evidence": "intercept only" if no_evidence
                    else (f"{n_fired} modelled mutation(s)" if n_fired else "not reported"),
        "basis": "probability margin, capped by the bundle's tier and by the "
                 "evidence the call rests on",
    }


def score(bundle: ModelBundle, features: dict[str, Any],
          threshold: float = 0.5) -> dict:
    """Score one isolate.

    Absent features default deliberately, not incidentally: a missing ``raw__``
    mutation feature is 0, meaning the variant is not carried, while a missing
    ``cov__`` covariate is mean-imputed to z = 0 so that a genotype-only call is
    not swamped by an unknown or out-of-distribution quality covariate.
    """
    m = bundle.model
    names, mean, scale, coef = (m["features"], m["scaler_mean"],
                                m["scaler_scale"], m["coef"])
    logit = m["intercept"]
    contributions = []
    for i, name in enumerate(names):
        if name in features:
            x = float(features[name])
        elif name.startswith("cov__"):
            x = mean[i]
        else:
            x = 0.0
        z = (x - mean[i]) / scale[i] if scale[i] else 0.0
        c = coef[i] * z
        logit += c
        if x != 0 or coef[i]:
            contributions.append({
                "feature": name, "value": x, "coef": round(coef[i], 4),
                "logit_contribution": round(c, 4),
                "direction": "R+" if coef[i] > 0 else "S-",
            })
    p = 1.0 / (1.0 + math.exp(-logit))
    contributions.sort(key=lambda r: -abs(r["logit_contribution"]))
    rng = bundle.operating_range
    return {
        "drug": bundle.drug,
        "prediction": "R" if p >= threshold else "S",
        "p_resistant": round(p, 4),
        "threshold": threshold,
        "operating_range": {
            "auc": rng.get("auc"), "tier": rng.get("tier"),
            "metric": rng.get("primary_metric"),
        },
        "reliability": reliability(
            p, rng.get("tier"),
            n_fired=sum(1 for c in contributions
                        if c["feature"].startswith("raw__") and c["value"] != 0)),
        "contributions": contributions,
    }
