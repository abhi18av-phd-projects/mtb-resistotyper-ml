"""Resistotyping service: VCF in, predictions out, two layers deep.

Layer 1 is the curated WHO catalogue read by piezo. Layer 2 is this project's
models. They do not vote: a drug the catalogue can grade is answered by the
catalogue and the model is not consulted for it at all, because a curated
assessment backed by phenotype evidence is not improved by averaging it with a
logistic regression. The models answer only where the catalogue is silent, which
is the role the paper claims for them.

Every row says which layer produced it. A Layer 2 row also carries the operating
range of the bundle that produced it, read from that bundle's card rather than
from a table in this service, so a model retrained to a different range is
reported with its own range and a bundle declaring none is reported as
unvalidated rather than silently assigned one.

If the catalogue cannot be loaded the service still answers, from models alone,
and says so: `catalogue_consulted` is false. That distinction matters. A model
call made after the catalogue declined to grade a variant is a different claim
from a model call made without asking.
"""

from __future__ import annotations

import csv
import io
import json
import os
import tempfile
import zipfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, Response

from mtb_resistotyper_ml import (__version__, ModelBundle, catalogue, explain,
                                 fired, score, vectorize)
from mtb_resistotyper_ml.report import render

from vcf_to_garc import ReferenceUnavailable, sample_ids, to_instance, variants_from_vcf

MODELS = Path(os.environ.get("MTB_MODELS", "/opt/models"))
REFERENCE = Path(os.environ.get("MTB_REFERENCE_GENBANK", "/opt/reference/NC_000962.3.gbk"))

app = FastAPI(title="mtb-resistotyper-ml", version="0.1.0")

DOWNLOAD_README = """mtb-resistotyper-ml predictions
===============================

predictions.tsv   one row per drug, with the reasoning flattened into a column
predictions.json  the full structured result, including every contribution
reasoning.txt     the same rendering `mtb-resistotyper-ml predict` prints

RESEARCH USE ONLY. NOT A DIAGNOSTIC. These models have not been prospectively
validated and have no regulatory clearance. Do not use any value in these files
to guide the treatment of a patient.

Each row names the layer that produced it. Layer 1 rows are graded by the WHO
catalogue (via piezo) and no model was consulted for that drug. Layer 2 rows are
this project's models, answering for variants the catalogue could not grade. The
`catalogue_consulted` column says whether Layer 1 actually ran; if it is false,
the catalogue was unavailable and the row rests on the model alone.

Each row carries the operating range declared by the bundle that produced it,
and a reliability band from 1 to 5. The band is the margin between the two class
probabilities, capped by that bundle's tier: a margin says how far apart this
model held the classes on this isolate, not whether the model is any good, so
the cap keeps a confident-looking margin from an indifferent model in its place.

Research use only. Not prospectively validated.

Reproduce any row with the same runner:
    pip install mtb-resistotyper-ml
    mtb-resistotyper-ml predict --input <instance.json> --model <bundle-dir>
"""


def available_drugs() -> list[str]:
    if not MODELS.exists():
        return []
    return sorted(d.name for d in MODELS.iterdir() if (d / "model.json").exists())


def score_instance(instance: dict) -> list[dict]:
    """Score one isolate against every deployed bundle.

    The result objects are byte-for-byte the shape ``mtb-resistotyper-ml predict
    --json`` produces, and each is rendered with the same ``report.render`` the
    command line prints. The service and the command line are then two front
    ends over one reasoning path rather than two implementations of it, so a
    laboratory cannot be shown one explanation on the web and another at a
    prompt for the same isolate.
    """
    # Layer 1 first. Everything the catalogue can grade is answered by the
    # catalogue, and the model is not consulted for that drug at all.
    variants = instance.get("variants", [])
    try:
        cat_available = True
        cat_version = catalogue.version()
    except Exception:
        cat_available = False
        cat_version = None

    results = []
    for drug in available_drugs():
        if cat_available:
            try:
                graded = catalogue.grade(variants, drug)
            except catalogue.CatalogueUnavailable:
                cat_available, graded = False, {"prediction": None}
            if graded["prediction"] is not None:
                results.append({
                    "sample_id": instance["sample_id"],
                    "lineage": instance.get("covariates", {}).get("lineage"),
                    "drug": drug,
                    "prediction": graded["prediction"],
                    "layer": 1,
                    "source": "WHO catalogue",
                    "n_mutations_carried": len(graded["graded"]),
                    "primary_drivers": [h["mutation"] for h in graded["graded"]
                                        if h["call"] == "R"],
                    "passengers_flagged": [],
                    "reasons": [{"mutation": h["mutation"], "call": h["call"],
                                 "role": "catalogue-graded", "on_target": True}
                                for h in graded["graded"]],
                    "ungraded_variants": graded["ungraded"],
                    "reliability": {"band": 5, "of": 5, "margin": None,
                                    "capped_by_tier": False, "capped_by_evidence": False,
                                    "evidence": "curated catalogue entry",
                                    "basis": "graded by the WHO catalogue, not modelled"},
                    "confidence": {"tier": "catalogue", "auc": None,
                                   "caveat": "Graded by "
                                   f"{graded['catalogue']['catalogue']} "
                                   f"{graded['catalogue']['version']}. No model was "
                                   "consulted for this drug."},
                    "provenance": {"runner": f"mtb-resistotyper-ml {__version__}",
                                   "catalogue": graded["catalogue"],
                                   "layer": 1, "catalogue_consulted": True},
                })
                continue
        bundle = ModelBundle.load(MODELS / drug)
        if not bundle.schema:
            continue
        row = vectorize(instance, bundle.schema)
        scored = score(bundle, row)
        results.append({
            "sample_id": instance["sample_id"],
            "lineage": instance.get("covariates", {}).get("lineage"),
            "n_mutations_carried": len(fired(row)),
            **scored,
            **explain(scored, bundle.drug),
            "layer": 2,
            "source": "model",
            "provenance": {
                "runner": f"mtb-resistotyper-ml {__version__}",
                "model_card": bundle.card.get("provenance", {}),
                "layer": 2,
                # True once Layer 1 ran and returned nothing gradable: the model
                # is then answering a question the catalogue declined, which is
                # a materially different claim from never having asked.
                "catalogue_consulted": cat_available,
                "catalogue": cat_version,
            },
        })
    return results


def rows_to_tsv(results: list[dict]) -> str:
    """One row per drug. The reasoning is flattened, not dropped: every carried
    mutation appears with its signed logit contribution, so the table can be read
    without the JSON beside it."""
    cols = ["sample_id", "drug", "prediction", "source", "layer", "p_resistant",
            "reliability_band", "reliability_margin", "tier", "auc",
            "n_mutations_carried", "primary_drivers", "passengers_flagged",
            "reasoning", "catalogue_consulted"]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in results:
        rel, conf = r["reliability"], r["confidence"]
        w.writerow({
            "sample_id": r["sample_id"], "drug": r["drug"],
            "prediction": r["prediction"], "p_resistant": r.get("p_resistant", ""),
            "reliability_band": f"{rel['band']}/{rel['of']}",
            "reliability_margin": rel["margin"] if rel["margin"] is not None else "",
            "tier": conf.get("tier") or "undeclared", "auc": conf.get("auc"),
            "n_mutations_carried": r["n_mutations_carried"],
            "primary_drivers": ";".join(r["primary_drivers"]) or "-",
            "passengers_flagged": ";".join(r["passengers_flagged"]) or "-",
            "reasoning": ";".join(
                f"{x['mutation']}({x['role']},logit{x['logit_contribution']:+.3f})"
                for x in r["reasons"]) or "intercept only",
            "layer": r.get("layer", 2),
            "source": r.get("source", "model"),
            "catalogue_consulted": r["provenance"].get("catalogue_consulted", False),
        })
    return buf.getvalue()


def rows_to_text(results: list[dict]) -> str:
    """The command line's own rendering, one block per drug."""
    return ("\n\n" + "-" * 72 + "\n\n").join(render(r) for r in results)


@app.get("/health")
def health() -> dict:
    try:
        cat = catalogue.version()
        catalogue.grade([], "RIF")          # forces the load
        cat_ok = True
    except Exception as exc:
        cat, cat_ok = {"error": str(exc)}, False
    return {"status": "ok", "drugs": available_drugs(),
            "reference_present": REFERENCE.exists(),
            "catalogue_available": cat_ok, "catalogue": cat}


@app.post("/predict")
async def predict(vcf: UploadFile | None = File(default=None),
                  instance_json: str | None = Form(default=None),
                  lineage: str = Form(default="unknown"),
                  fmt: str = Form(default="json")):
    if not available_drugs():
        raise HTTPException(503, f"no model bundles under {MODELS}")

    instances: list[dict] = []
    if instance_json:
        payload = json.loads(instance_json)
        instances = payload if isinstance(payload, list) else [payload]
    elif vcf is not None:
        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as fh:
            fh.write(await vcf.read())
            tmp = Path(fh.name)
        try:
            names = sample_ids(tmp) or [Path(vcf.filename or "sample").stem]
            variants = variants_from_vcf(tmp, REFERENCE)
            # One VCF may hold several samples. gumpy resolves one genome at a
            # time, so a multi-sample file is reported against its first column
            # and the rest are named but not scored, rather than silently
            # attributed the first sample's variants.
            instances = [to_instance(names[0], variants, lineage=lineage)]
            if len(names) > 1:
                instances[0]["_unscored_samples"] = names[1:]
        except ReferenceUnavailable as exc:
            raise HTTPException(503, str(exc)) from exc
        finally:
            tmp.unlink(missing_ok=True)
    else:
        raise HTTPException(400, "supply either a vcf file or instance_json")

    rows = [r for inst in instances for r in score_instance(inst)]

    if fmt == "txt":
        return Response(rows_to_text(rows), media_type="text/plain",
                        headers={"Content-Disposition": "attachment; filename=reasoning.txt"})
    if fmt == "tsv":
        return Response(rows_to_tsv(rows), media_type="text/tab-separated-values",
                        headers={"Content-Disposition": "attachment; filename=predictions.tsv"})
    if fmt == "zip":
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("predictions.tsv", rows_to_tsv(rows))
            z.writestr("predictions.json", json.dumps(rows, indent=2))
            z.writestr("reasoning.txt", rows_to_text(rows))
            z.writestr("README.txt", DOWNLOAD_README)
        return Response(buf.getvalue(), media_type="application/zip",
                        headers={"Content-Disposition": "attachment; filename=predictions.zip"})
    return JSONResponse(rows)


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (Path(__file__).with_name("static") / "index.html").read_text()
