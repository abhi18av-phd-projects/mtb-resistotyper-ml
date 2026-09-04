"""Resistotyping service: a VCF the cluster already holds, in; predictions out.

Two layers deep. Layer 1 is the curated WHO catalogue read by piezo; Layer 2 is
this project's models. They do not vote: a drug the catalogue can grade is
answered by the catalogue and the model is not consulted for it, because a
curated assessment backed by phenotype evidence is not improved by averaging it
with a logistic regression. Measured through piezo, 146 of the 223 mutations the
deployed bundles model lie outside the catalogue's grading, so Layer 2's remit
is most of what the models describe rather than an edge case.

The primary input is an object the researcher already uploaded with `abc data
upload` or through the browser upload path, not a second copy posted to this
service. Direct upload remains as a fallback for someone with a file and no
cluster account.

Both front ends -- this service and the `mtb-resistotyper-ml` command -- call
`mtb_resistotyper_ml.layers.resolve_all` and render with the same
`report.render`, so neither can drift into telling a different story about the
same isolate.
"""

from __future__ import annotations

import asyncio
import contextlib
import csv
import io
import json
import os
import tempfile
import zipfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, Response

from mtb_resistotyper_ml import __version__, catalogue, catalogue_ready, resolve_all
from mtb_resistotyper_ml.report import render

import normalise
import storage
from report_html import render_report
from vcf_to_garc import (ReferenceUnavailable, genes_of_interest, to_instance,
                         variants_from_vcf)

MODELS = Path(os.environ.get("MTB_MODELS", "/opt/models"))
REFERENCE = Path(os.environ.get("MTB_REFERENCE_GENBANK", "/opt/reference/NC_000962.3.gbk"))
HERE = Path(__file__).parent

app = FastAPI(title="mtb-resistotyper-ml", version=__version__)

_SWEEP_STATE: dict = {"last": None, "result": None, "lifecycle": None}


async def _retention_loop() -> None:
    """Apply the expiry rule, then sweep daily.

    The rule is primary: it is server-side, so it keeps expiring uploads even if
    this process is stopped, redeployed, or crash-looping. The sweep is the
    backstop, because a lifecycle rule that was never applied fails silently and
    looks exactly like one that is working. Both report through /health, so the
    difference is visible without reading a log.
    """
    while True:
        try:
            _SWEEP_STATE["lifecycle"] = storage.ensure_lifecycle()
        except Exception as exc:
            _SWEEP_STATE["lifecycle"] = {"error": str(exc)[:200]}
        try:
            _SWEEP_STATE["result"] = storage.sweep_anonymous()
        except Exception as exc:
            _SWEEP_STATE["result"] = {"error": str(exc)[:200]}
        _SWEEP_STATE["last"] = __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat()
        await asyncio.sleep(24 * 3600)


@app.on_event("startup")
async def _start_retention() -> None:
    if storage.anon_enabled():
        app.state.retention = asyncio.create_task(_retention_loop())


@app.on_event("shutdown")
async def _stop_retention() -> None:
    task = getattr(app.state, "retention", None)
    if task:
        task.cancel()
        with contextlib.suppress(asyncio.CancelledError):
            await task


@app.post("/admin/sweep")
def sweep(dry_run: bool = True) -> dict:
    """Run the retention sweep now. Defaults to a dry run so a probe cannot delete."""
    return storage.sweep_anonymous(dry_run=dry_run)

DOWNLOAD_README = """mtb-resistotyper-ml predictions
===============================

report.html         the whole result set as one self-contained page: the
                    susceptibility grid, per-drug reasoning, the contribution
                    plots and the provenance. Open this first.
predictions.tsv     one row per drug, with the reasoning flattened into a column
predictions.json    the full structured result, including every contribution
reasoning.txt       the same rendering `mtb-resistotyper-ml predict` prints
causal_chains.json  for each model-derived call, the chain from observed variant
                    to reported call, separating the contributors with an
                    established mechanism for the drug from those carrying weight
                    only through co-occurrence in the training population
inputs.json         which object was read and which variants were derived from it

RESEARCH USE ONLY. NOT A DIAGNOSTIC. The Layer 2 models have not been
prospectively validated and have no regulatory clearance. Do not use any value
in these files to guide the treatment of a patient.

Each row names the layer that produced it. Layer 1 rows are graded by the WHO
catalogue (via piezo) and no model was consulted for that drug. Layer 2 rows are
this project's models, answering for variants the catalogue could not grade. The
catalogue_consulted column says whether Layer 1 actually ran; if it is false the
catalogue was unavailable and the row rests on the model alone.

Layer 2 rows carry a reliability band from 1 to 5: the margin between the two
class probabilities, capped by the bundle's tier and by whether any modelled
mutation actually fired.

Reproduce any row with the same runner and the same inputs:
    pip install mtb-resistotyper-ml
    mtb-resistotyper-ml predict --input <instance.json> --models <bundles> \\
        --catalogue <who.csv> --outdir .
"""

TSV_COLUMNS = ["sample_id", "drug", "prediction", "source", "layer", "p_resistant",
               "reliability_band", "reliability_margin", "tier", "auc",
               "n_mutations_carried", "primary_drivers", "passengers_flagged",
               "reasoning", "catalogue_consulted", "catalogue_version"]


def available_drugs() -> list[str]:
    if not MODELS.exists():
        return []
    return sorted(d.name for d in MODELS.iterdir() if (d / "model.json").exists())


def _reasoning(r: dict) -> str:
    if not r.get("reasons"):
        return "catalogue: no graded mutation" if r["layer"] == 1 else "intercept only"
    if r["layer"] == 1:
        return ";".join(f"{x['mutation']}={x['call']}(WHO)" for x in r["reasons"])
    return ";".join(f"{x['mutation']}({x['role']},logit{x['logit_contribution']:+.3f})"
                    for x in r["reasons"])


def rows_to_tsv(results: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in results:
        rel, conf, prov = r["reliability"], r["confidence"], r["provenance"]
        cat = prov.get("catalogue") or {}
        w.writerow({
            "sample_id": r["sample_id"], "drug": r["drug"],
            "prediction": r["prediction"], "source": r["source"], "layer": r["layer"],
            "p_resistant": "" if r.get("p_resistant") is None else r["p_resistant"],
            "reliability_band": f"{rel['band']}/{rel['of']}",
            "reliability_margin": "" if rel["margin"] is None else rel["margin"],
            "tier": conf.get("tier") or "undeclared",
            "auc": "" if conf.get("auc") is None else conf["auc"],
            "n_mutations_carried": r["n_mutations_carried"],
            "primary_drivers": ";".join(r["primary_drivers"]) or "-",
            "passengers_flagged": ";".join(r.get("passengers_flagged", [])) or "-",
            "reasoning": _reasoning(r),
            "catalogue_consulted": prov.get("catalogue_consulted", False),
            "catalogue_version": cat.get("version") or "-",
        })
    return buf.getvalue()


def rows_to_text(results: list[dict]) -> str:
    return ("\n\n" + "-" * 72 + "\n\n").join(render(r) for r in results)


def _instances_from_vcf(path: Path, lineage: str, workdir: Path) -> list[dict]:
    """One ideal-input instance per sample in the file.

    A multi-sample cohort produces several. Earlier this service scored the first
    column and merely named the others, which reports one isolate's genotype
    under a shared filename; bcftools splits them properly instead.
    """
    genes = genes_of_interest(MODELS)
    out = []
    for sample, ready in normalise.normalise(path, workdir, genes=genes, genbank=REFERENCE):
        variants = variants_from_vcf(ready, REFERENCE, genes=genes)
        out.append(to_instance(sample, variants, lineage=lineage))
    return out


@app.get("/health")
def health() -> dict:
    try:
        cat_ok = catalogue_ready()
        cat = catalogue.version() if cat_ok else {"error": "unavailable"}
    except Exception as exc:
        cat_ok, cat = False, {"error": str(exc)}
    try:
        storage.list_vcfs(limit=1)
        store_ok = True
    except storage.ListingDisabled:
        store_ok = True          # listing is off by policy; storage itself is fine
    except Exception:
        store_ok = False
    return {"status": "ok", "version": __version__, "drugs": available_drugs(),
            "reference_present": REFERENCE.exists(),
            "bcftools": normalise.available(),
            "catalogue_available": cat_ok, "catalogue": cat,
            "storage_available": store_ok, "bucket": storage.BUCKET,
            "listing_enabled": bool(storage.LIST_PREFIXES),
            "list_prefixes": storage.LIST_PREFIXES,
            "anonymous_uploads": storage.anon_enabled(),
            "retention": {
                "hours": storage.ANON_RETENTION_HOURS,
                "prefix": storage.ANON_PREFIX,
                "max_upload_bytes": storage.ANON_MAX_BYTES,
                "bucket_rule": storage.lifecycle_state() if storage.anon_enabled() else None,
                "last_sweep": _SWEEP_STATE["last"],
                "last_sweep_result": _SWEEP_STATE["result"],
            }}


@app.get("/objects")
def objects() -> JSONResponse:
    """VCFs already uploaded to the group bucket."""
    try:
        return JSONResponse(storage.list_vcfs())
    except storage.ListingDisabled as exc:
        raise HTTPException(403, str(exc)) from exc
    except storage.StorageUnavailable as exc:
        raise HTTPException(503, str(exc)) from exc


@app.post("/predict")
async def predict(object_key: str | None = Form(default=None),
                  vcf: UploadFile | None = File(default=None),
                  instance_json: str | None = Form(default=None),
                  lineage: str = Form(default="unknown"),
                  fmt: str = Form(default="json")):
    if not available_drugs():
        raise HTTPException(503, f"no model bundles under {MODELS}")

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        source: dict = {}
        try:
            if instance_json:
                payload = json.loads(instance_json)
                instances = payload if isinstance(payload, list) else [payload]
                source = {"kind": "instance_json"}
            elif object_key:
                local = storage.fetch(object_key, work / "in")
                instances = _instances_from_vcf(local, lineage, work / "norm")
                source = {"kind": "object", "bucket": storage.BUCKET, "key": object_key}
            elif vcf is not None:
                raw = await vcf.read()
                local = work / "in" / (vcf.filename or "upload.vcf")
                local.parent.mkdir(parents=True, exist_ok=True)
                local.write_bytes(raw)
                source = {"kind": "upload", "filename": vcf.filename}
                if storage.anon_enabled():
                    # Retained briefly so a tester can re-run without re-uploading,
                    # and expired by the bucket rule so "briefly" is enforced by the
                    # store rather than promised by this code.
                    try:
                        key = storage.store_anonymous(raw, vcf.filename or "upload.vcf")
                        source |= {"stored_key": key,
                                   "retention_hours": storage.ANON_RETENTION_HOURS}
                    except storage.StorageUnavailable as exc:
                        source["not_stored"] = str(exc)
                instances = _instances_from_vcf(local, lineage, work / "norm")
            else:
                raise HTTPException(400, "supply object_key, a vcf file, or instance_json")
        except storage.StorageUnavailable as exc:
            raise HTTPException(503, f"storage: {exc}") from exc
        except normalise.NormalisationFailed as exc:
            raise HTTPException(422, f"VCF could not be normalised: {exc}") from exc
        except ReferenceUnavailable as exc:
            raise HTTPException(503, str(exc)) from exc

        rows = [r for inst in instances
                for r in resolve_all(inst, MODELS, use_catalogue=catalogue_ready())]

    if fmt == "html":
        return Response(render_report(rows, inputs={"source": source, "lineage": lineage}),
                        media_type="text/html")
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
            chains = {f"{r['sample_id']}/{r['drug']}": r["causal_chain"]
                      for r in rows if r.get("causal_chain")}
            if chains:
                z.writestr("causal_chains.json", json.dumps(chains, indent=2))
            z.writestr("inputs.json", json.dumps(
                {"source": source, "lineage": lineage, "runner": __version__,
                 "samples": [{"sample_id": i["sample_id"],
                              "variants": i.get("variants", [])} for i in instances]},
                indent=2, default=str))
            z.writestr("README.txt", DOWNLOAD_README)
            # The report is the artefact a reader actually opens; the rest are
            # what they reach for once it has told them where to look.
            z.writestr("report.html", render_report(
                rows, inputs={"source": source, "lineage": lineage}))
        return Response(buf.getvalue(), media_type="application/zip",
                        headers={"Content-Disposition": "attachment; filename=predictions.zip"})
    return JSONResponse(rows)


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (HERE / "static" / "index.html").read_text()
