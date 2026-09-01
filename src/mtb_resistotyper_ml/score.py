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
        "contributions": contributions,
    }
