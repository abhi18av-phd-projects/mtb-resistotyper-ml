"""Resistotyping service: VCF in, predictions out.

Layer 2 only. These models are meant to fire where a curated catalogue returns
Unknown or Fail, and the catalogue layer is not wired in here, so every call this
service makes is a model call and is labelled as such. Nothing it returns
overrides a catalogue, because it never consults one.

Every prediction carries the operating range of the bundle that produced it,
read from that bundle's card rather than from a table in this service. A model
retrained to a different range is then reported with its own range, and a bundle
declaring none is reported as unvalidated rather than silently assigned one.
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

from mtb_resistotyper_ml import __version__, ModelBundle, explain, fired, score, vectorize
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

These are Layer 2 model calls. No curated WHO catalogue was consulted, so
nothing here overrides a catalogue result; the models are intended for variants
a catalogue grades Unknown or Fail.

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
    results = []
    for drug in available_drugs():
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
            "provenance": {
                "runner": f"mtb-resistotyper-ml {__version__}",
                "model_card": bundle.card.get("provenance", {}),
                "layer": 2,
                "catalogue_consulted": False,
            },
        })
    return results


def rows_to_tsv(results: list[dict]) -> str:
    """One row per drug. The reasoning is flattened, not dropped: every carried
    mutation appears with its signed logit contribution, so the table can be read
    without the JSON beside it."""
    cols = ["sample_id", "drug", "prediction", "p_resistant", "reliability_band",
            "reliability_margin", "tier", "auc", "n_mutations_carried",
            "primary_drivers", "passengers_flagged", "reasoning",
            "layer", "catalogue_consulted"]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in results:
        rel, conf = r["reliability"], r["confidence"]
        w.writerow({
            "sample_id": r["sample_id"], "drug": r["drug"],
            "prediction": r["prediction"], "p_resistant": r["p_resistant"],
            "reliability_band": f"{rel['band']}/{rel['of']}",
            "reliability_margin": rel["margin"],
            "tier": conf.get("tier") or "undeclared", "auc": conf.get("auc"),
            "n_mutations_carried": r["n_mutations_carried"],
            "primary_drivers": ";".join(r["primary_drivers"]) or "-",
            "passengers_flagged": ";".join(r["passengers_flagged"]) or "-",
            "reasoning": ";".join(
                f"{x['mutation']}({x['role']},logit{x['logit_contribution']:+.3f})"
                for x in r["reasons"]) or "intercept only",
            "layer": 2, "catalogue_consulted": False,
        })
    return buf.getvalue()


def rows_to_text(results: list[dict]) -> str:
    """The command line's own rendering, one block per drug."""
    return ("\n\n" + "-" * 72 + "\n\n").join(render(r) for r in results)


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "drugs": available_drugs(),
            "reference_present": REFERENCE.exists()}


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
