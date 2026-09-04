"""The layer boundary, in one place.

The command line and the web service must not disagree about where an answer
came from. Putting the boundary here rather than in either front end is what
makes that structural: both call ``resolve`` and neither decides for itself when
a model may fire.

Layer 1 is the curated WHO catalogue read by piezo. Layer 2 is this project's
models. They do not vote. A drug the catalogue grades is answered by the
catalogue and the model is not consulted for it, because a curated assessment
backed by phenotype evidence is not improved by averaging it with a logistic
regression.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from mtb_resistotyper_ml import catalogue as _catalogue
from mtb_resistotyper_ml.explain import causal_chain, explain
from mtb_resistotyper_ml.score import ModelBundle, score
from mtb_resistotyper_ml.vectorize import fired, vectorize
from mtb_resistotyper_ml.version import __version__


def catalogue_ready() -> bool:
    try:
        _catalogue.grade([], "RIF")
        return True
    except Exception:
        return False


def resolve(instance: dict, bundle: ModelBundle, *, threshold: float = 0.5,
            use_catalogue: bool = True) -> dict:
    """Answer one drug for one isolate, from whichever layer owns it."""
    drug = bundle.drug
    variants = instance.get("variants", [])
    common: dict[str, Any] = {
        "sample_id": instance.get("sample_id"),
        "lineage": instance.get("covariates", {}).get("lineage"),
        "drug": drug,
    }

    if use_catalogue:
        try:
            graded = _catalogue.grade(variants, drug)
        except _catalogue.CatalogueUnavailable:
            use_catalogue, graded = False, {"prediction": None}
        if graded["prediction"] is not None:
            return {
                **common,
                "prediction": graded["prediction"],
                "p_resistant": None,
                "layer": 1,
                "source": "WHO catalogue",
                "n_mutations_carried": len(graded["graded"]),
                "primary_drivers": [h["mutation"] for h in graded["graded"]
                                    if h["call"] == "R"],
                "passengers_flagged": [],
                "reasons": [{"mutation": h["mutation"], "call": h["call"],
                             "role": "catalogue-graded", "on_target": True,
                             "coefficient": None, "logit_contribution": None}
                            for h in graded["graded"]],
                "ungraded_variants": graded["ungraded"],
                "reliability": {"band": 5, "of": 5, "margin": None,
                                "capped_by_tier": False, "capped_by_evidence": False,
                                "evidence": "curated catalogue entry",
                                "basis": "graded by the WHO catalogue, not modelled"},
                "confidence": {"tier": "catalogue", "auc": None,
                               "caveat": f"Graded by {graded['catalogue']['catalogue']} "
                                         f"{graded['catalogue']['version']}. No model "
                                         f"was consulted for this drug."},
                "provenance": {"runner": f"mtb-resistotyper-ml {__version__}",
                               "catalogue": graded["catalogue"],
                               "layer": 1, "catalogue_consulted": True},
            }

    if not bundle.schema:
        raise ValueError(f"{drug}: no feature_schema.json; the input contract is unknown")
    row = vectorize(instance, bundle.schema)
    scored = score(bundle, row, threshold=threshold)
    result = {
        **common,
        "n_mutations_carried": len(fired(row)),
        **scored,
        **explain(scored, drug),
        "layer": 2,
        "source": "model",
        "provenance": {
            "runner": f"mtb-resistotyper-ml {__version__}",
            "model_card": bundle.card.get("provenance", {}),
            "layer": 2,
            # True once Layer 1 ran and declined: a model answering a question the
            # catalogue could not grade is a different claim from one answering
            # without asking.
            "catalogue_consulted": use_catalogue,
            "catalogue": _catalogue.version() if use_catalogue else None,
        },
    }
    # Only model rows carry a chain. A catalogue grade is not a chain of
    # inference from this tool's premises; it is a citation of somebody else's
    # curated judgement, and dressing it up as reasoning would misrepresent both.
    result["causal_chain"] = causal_chain(result, bundle.card)
    return result


def resolve_all(instance: dict, models_dir: Path, *, threshold: float = 0.5,
                use_catalogue: bool = True, drugs: list[str] | None = None) -> list[dict]:
    """Answer every deployed drug for one isolate."""
    available = sorted(d.name for d in Path(models_dir).iterdir()
                       if (d / "model.json").exists())
    out = []
    for drug in (drugs or available):
        d = Path(models_dir) / drug
        if not (d / "model.json").exists():
            continue
        out.append(resolve(instance, ModelBundle.load(d), threshold=threshold,
                           use_catalogue=use_catalogue))
    return out
