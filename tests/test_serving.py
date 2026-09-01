"""Golden tests for the serving path.

These are the tests a reviewer runs to confirm the tool does what the paper says.
Each fixture in ``tests/data`` has a recorded expected output in ``tests/expected``;
the point of committing both is that a change in scoring cannot pass silently.

Set MTB_MODELS to a checkout of mtb-resistotyper-ml-models to run them.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from mtb_resistotyper_ml import ModelBundle, explain, fired, score, vectorize
from mtb_resistotyper_ml.report import render

HERE = Path(__file__).parent
MODELS = Path(os.environ.get("MTB_MODELS", HERE.parent.parent / "mtb-resistotyper-ml-models"))
RIF = MODELS / "models" / "RIF"
needs_models = pytest.mark.skipif(not RIF.exists(), reason="set MTB_MODELS to a models checkout")

CASES = ["rif-resistant", "rif-susceptible"]


def _run(case: str) -> dict:
    bundle = ModelBundle.load(RIF)
    instance = json.loads((HERE / "data" / f"{case}.input.json").read_text())
    row = vectorize(instance, bundle.schema)
    scored = score(bundle, row)
    return {"sample_id": instance["sample_id"],
            "lineage": instance["covariates"]["lineage"],
            "n_mutations_carried": len(fired(row)),
            **scored, **explain(scored, bundle.drug)}


@needs_models
@pytest.mark.parametrize("case", CASES)
def test_matches_expected_output(case: str) -> None:
    got = _run(case)
    want = json.loads((HERE / "expected" / f"{case}.json").read_text())
    for key in ("prediction", "p_resistant", "n_mutations_carried", "primary_drivers"):
        assert got[key] == want[key], key


@needs_models
def test_operating_range_travels_with_the_model() -> None:
    """The tier must come from the bundle, never from a table inside this package.

    A hardcoded tier goes stale the moment a model is retrained, and would then
    misreport the range of an artefact that correctly declares its own.
    """
    bundle = ModelBundle.load(RIF)
    rng = bundle.operating_range
    assert rng["tier"] in {"usable", "moderate"}
    assert 0.5 < rng["auc"] <= 1.0
    assert "lineage" in rng["primary_metric"]

    bundle.card = {}
    out = explain(score(bundle, {}), bundle.drug)
    assert out["confidence"]["tier"] is None
    assert "unvalidated" in out["confidence"]["caveat"]


@needs_models
def test_contributions_sum_to_the_logit() -> None:
    """The explanation IS the model, so the parts must add up to the whole."""
    import math
    bundle = ModelBundle.load(RIF)
    instance = json.loads((HERE / "data" / "rif-resistant.input.json").read_text())
    scored = score(bundle, vectorize(instance, bundle.schema))
    total = bundle.model["intercept"] + sum(c["logit_contribution"]
                                            for c in scored["contributions"])
    assert scored["p_resistant"] == pytest.approx(1 / (1 + math.exp(-total)), abs=1e-3)


@needs_models
def test_absent_variant_is_wild_type_not_missing() -> None:
    bundle = ModelBundle.load(RIF)
    empty = vectorize({"variants": [], "covariates": {}}, bundle.schema)
    assert all(v == 0 for k, v in empty.items() if k.startswith("raw__"))
    assert score(bundle, empty)["prediction"] == "S"


@needs_models
@pytest.mark.parametrize("case", CASES)
def test_report_renders(case: str) -> None:
    text = render(_run(case))
    assert (HERE / "expected" / f"{case}.txt").read_text().strip() == text.strip()
