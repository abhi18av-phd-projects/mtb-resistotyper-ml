"""Command-line entry point.

Shaped for a workflow engine as much as for a person. A Nextflow process wants
one command that takes a file in, writes files out with predictable names, exits
non-zero on failure, and never prompts: hence ``--outdir`` with fixed basenames,
``--tsv`` for a channel-friendly table, and ``--quiet``. It also wants to pin what
it ran, so ``--catalogue`` and ``--models`` are explicit paths rather than
discovered state, and every output records both.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

from mtb_resistotyper_ml import catalogue as _catalogue
from mtb_resistotyper_ml.ingest import BcftoolsUnavailable, ReferenceUnavailable, instances_from_vcf
from mtb_resistotyper_ml.layers import catalogue_ready, resolve_all
from mtb_resistotyper_ml.report import render
from mtb_resistotyper_ml.score import ModelBundle
from mtb_resistotyper_ml.version import __version__

TSV_COLUMNS = ["sample_id", "drug", "prediction", "source", "layer", "p_resistant",
               "reliability_band", "reliability_margin", "tier", "auc",
               "n_mutations_carried", "primary_drivers", "passengers_flagged",
               "reasoning", "catalogue_consulted", "catalogue_version", "runner"]


def _rows_to_tsv(results: list[dict], fh) -> None:
    w = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
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
            "runner": prov.get("runner", ""),
        })


def _reasoning(r: dict) -> str:
    """Why this row says what it says, flattened but not lost.

    A table that gives a call without its grounds forces the reader back to the
    JSON, and in practice they do not go. Layer 1 rows carry the graded mutations
    and their catalogue calls; Layer 2 rows carry each mutation's signed logit
    contribution, which for this model class IS the explanation rather than an
    approximation of it.
    """
    if not r.get("reasons"):
        return "catalogue: no graded mutation" if r["layer"] == 1 else "intercept only"
    if r["layer"] == 1:
        return ";".join(f"{x['mutation']}={x['call']}(WHO)" for x in r["reasons"])
    return ";".join(f"{x['mutation']}({x['role']},logit{x['logit_contribution']:+.3f})"
                    for x in r["reasons"])


def _load_instance(path: str) -> dict:
    return json.loads(Path(path).read_text())


def _predict(args: argparse.Namespace) -> int:
    if args.catalogue:
        _catalogue.CATALOGUE = Path(args.catalogue)
    use_cat = not args.no_catalogue
    if use_cat and not catalogue_ready():
        # Loud, not silent. A run that quietly skipped Layer 1 produces model
        # calls for drugs the catalogue would have graded, which is a different
        # analysis wearing the same output shape.
        msg = (f"Layer 1 unavailable (piezo missing, or no catalogue at "
               f"{_catalogue.CATALOGUE}). ")
        if args.require_catalogue:
            sys.exit(msg + "Refusing: --require-catalogue was given.")
        print(msg + "Continuing with models only; rows will say "
                    "catalogue_consulted=false.", file=sys.stderr)

    models = Path(args.models)
    if args.vcf:
        try:
            instances = instances_from_vcf(args.vcf, args.reference, models=models,
                                           lineage=args.lineage or "unknown")
        except (ReferenceUnavailable, BcftoolsUnavailable) as exc:
            sys.exit(str(exc))
        if not instances:
            sys.exit(f"no samples found in {args.vcf}")
    else:
        instances = [_load_instance(args.input)]
        if args.lineage:
            instances[0].setdefault("covariates", {})["lineage"] = args.lineage
    if args.model:                     # single-bundle mode, as before
        bundles = [Path(args.model)]
        models, drugs = bundles[0].parent, [bundles[0].name]
    else:
        drugs = args.drugs.split(",") if args.drugs else None

    every: list[dict] = []
    for instance in instances:
        results = resolve_all(instance, models, threshold=args.threshold,
                              use_catalogue=use_cat, drugs=drugs)
        if not results:
            sys.exit(f"no model bundles under {models}")
        every.extend(results)

        if args.outdir:
            # One set of files per sample. A cohort file scored into a single
            # merged output would need the reader to re-separate the isolates,
            # which is exactly the confusion the per-sample split exists to avoid.
            out = Path(args.outdir)
            out.mkdir(parents=True, exist_ok=True)
            sid = instance.get("sample_id", "sample")
            (out / f"{sid}.predictions.json").write_text(json.dumps(results, indent=2) + "\n")
            with (out / f"{sid}.predictions.tsv").open("w") as fh:
                _rows_to_tsv(results, fh)
            (out / f"{sid}.reasoning.txt").write_text(
                ("\n\n" + "-" * 72 + "\n\n").join(render(r) for r in results) + "\n")
            if not args.quiet:
                print(f"wrote {out}/{sid}.predictions.{{json,tsv}} and {sid}.reasoning.txt")

    if args.outdir:
        return 0
    if args.json:
        print(json.dumps(every, indent=2))
    elif args.tsv:
        _rows_to_tsv(every, sys.stdout)
    else:
        print(("\n\n" + "-" * 72 + "\n\n").join(render(r) for r in every))
    return 0


def _describe(args: argparse.Namespace) -> int:
    """Print what a model bundle declares about itself, without scoring anything."""
    bundle = ModelBundle.load(args.model)
    rng = bundle.operating_range
    print(f"drug        {bundle.drug}")
    print(f"model       {bundle.model.get('model')}")
    print(f"features    {len(bundle.model['features'])}")
    if rng:
        print(f"AUC         {rng.get('auc')}   ({rng.get('primary_metric')})")
        print(f"tier        {rng.get('tier')}")
        per = rng.get("per_lineage_auc") or {}
        if per:
            print("per lineage " + "  ".join(f"{k}={v}" for k, v in sorted(per.items())))
        print(f"caveat      {rng.get('caveat', '')}")
    else:
        print("AUC         UNDECLARED — this bundle states no operating range")
    return 0


def _layers(args: argparse.Namespace) -> int:
    """Report which layer would answer each drug, without scoring an isolate.

    Written for provenance rather than curiosity: a pipeline that records this
    alongside its results can say afterwards which calls a catalogue edition was
    responsible for, without re-running anything.
    """
    if args.catalogue:
        _catalogue.CATALOGUE = Path(args.catalogue)
    ready = catalogue_ready()
    print(f"catalogue   {'available' if ready else 'UNAVAILABLE'}")
    if ready:
        for k, v in _catalogue.version().items():
            print(f"  {k:9} {v}")
    models = Path(args.models)
    print(f"models      {models}")
    for d in sorted(p for p in models.iterdir() if (p / "model.json").exists()):
        b = ModelBundle.load(d)
        n = len([f for f in b.schema.get("features", []) if f["kind"] == "mutation"])
        print(f"  {b.drug:5} {n:3d} modelled mutations   tier={b.operating_range.get('tier')}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mtb-resistotyper-ml",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Two-layer M. tuberculosis resistance calls from variant calls.\n"
                    "Layer 1 is the curated WHO catalogue (piezo); Layer 2 is a versioned\n"
                    "model bundle. They do not vote: a drug the catalogue grades is answered\n"
                    "by the catalogue and no model is consulted for it.",
        epilog="Nextflow:\n"
               "  mtb-resistotyper-ml predict --input ${instance} --models ${models} \\\n"
               "      --catalogue ${cat} --outdir . --quiet\n"
               "  // emits <sample>.predictions.{json,tsv} and <sample>.reasoning.txt\n")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("predict", help="score one isolate across every deployed drug")
    src = pr.add_mutually_exclusive_group(required=True)
    src.add_argument("--vcf", help="VCF/gVCF/BCF, single- or multi-sample (e.g. MAGMA output). "
                                   "Normalised with bcftools and converted to GARC, the same "
                                   "path the web service takes")
    src.add_argument("--input", help="ideal-input JSON (GARC variants + covariates), "
                                     "for an isolate already converted")
    pr.add_argument("--reference", default=os.environ.get(
                        "MTB_REFERENCE_GENBANK", "/opt/reference/NC_000962.3.gbk"),
                    help="H37Rv GenBank for --vcf (env: MTB_REFERENCE_GENBANK)")
    pr.add_argument("--models", default=os.environ.get("MTB_MODELS", "models"),
                    help="directory of model bundles (env: MTB_MODELS)")
    pr.add_argument("--model", help="score a single bundle directory instead of all of them")
    pr.add_argument("--drugs", help="comma-separated subset, e.g. RIF,INH")
    pr.add_argument("--catalogue", default=os.environ.get("MTB_CATALOGUE"),
                    help="WHO GARC1 catalogue CSV for Layer 1 (env: MTB_CATALOGUE)")
    pr.add_argument("--no-catalogue", action="store_true",
                    help="skip Layer 1 entirely; every call becomes a model call")
    pr.add_argument("--require-catalogue", action="store_true",
                    help="exit non-zero if Layer 1 cannot run, instead of continuing")
    pr.add_argument("--lineage", help="override the instance's lineage covariate")
    pr.add_argument("--threshold", type=float, default=0.5)
    pr.add_argument("--outdir", help="write <sample>.predictions.{json,tsv} + .reasoning.txt here")
    pr.add_argument("--json", action="store_true", help="structured result to stdout")
    pr.add_argument("--tsv", action="store_true", help="one row per drug to stdout")
    pr.add_argument("--quiet", action="store_true", help="suppress progress chatter")
    pr.set_defaults(func=_predict)

    de = sub.add_parser("describe", help="print a bundle's declared operating range")
    de.add_argument("--model", required=True, help="model bundle directory")
    de.set_defaults(func=_describe)

    la = sub.add_parser("layers", help="report which layers are available and what they cover")
    la.add_argument("--models", default=os.environ.get("MTB_MODELS", "models"))
    la.add_argument("--catalogue", default=os.environ.get("MTB_CATALOGUE"))
    la.set_defaults(func=_layers)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
