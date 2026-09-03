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

from mtb_resistotyper_ml import ModelBundle, explain, fired, score, vectorize

from vcf_to_garc import ReferenceUnavailable, sample_ids, to_instance, variants_from_vcf

MODELS = Path(os.environ.get("MTB_MODELS", "/opt/models"))
REFERENCE = Path(os.environ.get("MTB_REFERENCE_GENBANK", "/opt/reference/NC_000962.3.gbk"))

app = FastAPI(title="mtb-resistotyper-ml", version="0.1.0")


def available_drugs() -> list[str]:
    if not MODELS.exists():
        return []
    return sorted(d.name for d in MODELS.iterdir() if (d / "model.json").exists())


def score_instance(instance: dict) -> list[dict]:
    rows = []
    for drug in available_drugs():
        bundle = ModelBundle.load(MODELS / drug)
        if not bundle.schema:
            continue
        row_features = vectorize(instance, bundle.schema)
        scored = score(bundle, row_features)
        reasoning = explain(scored, drug)
        rng = scored.get("operating_range", {})
        rows.append({
            "sample_id": instance["sample_id"],
            "drug": drug,
            "prediction": scored["prediction"],
            "p_resistant": scored["p_resistant"],
            "tier": rng.get("tier"),
            "auc": rng.get("auc"),
            "metric": rng.get("metric"),
            "n_mutations_carried": len(fired(row_features)),
            "primary_drivers": reasoning["primary_drivers"],
            "passengers_flagged": reasoning["passengers_flagged"],
            "caveat": reasoning["confidence"]["caveat"],
            "layer": 2,
            "catalogue_consulted": False,
        })
    return rows


def rows_to_tsv(rows: list[dict]) -> str:
    cols = ["sample_id", "drug", "prediction", "p_resistant", "tier", "auc",
            "n_mutations_carried", "primary_drivers", "passengers_flagged",
            "layer", "catalogue_consulted"]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({**r,
                    "primary_drivers": ";".join(r["primary_drivers"]),
                    "passengers_flagged": ";".join(r["passengers_flagged"])})
    return buf.getvalue()


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

    if fmt == "tsv":
        return Response(rows_to_tsv(rows), media_type="text/tab-separated-values",
                        headers={"Content-Disposition": "attachment; filename=predictions.tsv"})
    if fmt == "zip":
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("predictions.tsv", rows_to_tsv(rows))
            z.writestr("predictions.json", json.dumps(rows, indent=2))
            z.writestr("README.txt",
                       "Layer 2 model predictions. No curated catalogue was consulted.\n"
                       "Each row carries the operating range of the bundle that produced it.\n"
                       "Research use only; not prospectively validated.\n")
        return Response(buf.getvalue(), media_type="application/zip",
                        headers={"Content-Disposition": "attachment; filename=predictions.zip"})
    return JSONResponse(rows)


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (Path(__file__).with_name("static") / "index.html").read_text()
