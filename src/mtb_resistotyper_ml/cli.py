"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from mtb_resistotyper_ml import __version__
from mtb_resistotyper_ml.explain import explain
from mtb_resistotyper_ml.report import render
from mtb_resistotyper_ml.score import ModelBundle, score
from mtb_resistotyper_ml.vectorize import fired, vectorize


def _predict(args: argparse.Namespace) -> int:
    bundle = ModelBundle.load(args.model)
    if not bundle.schema:
        sys.exit(f"{args.model}: no feature_schema.json; the input contract is unknown")
    instance = json.loads(Path(args.input).read_text())

    row = vectorize(instance, bundle.schema)
    scored = score(bundle, row, threshold=args.threshold)
    result = {
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
    }
    if args.json:
        out = json.dumps(result, indent=2)
    else:
        out = render(result)
    if args.out:
        Path(args.out).write_text(out + "\n")
    else:
        print(out)
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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mtb-resistotyper-ml",
        description="Predict M. tuberculosis drug resistance from variant calls, "
                    "using a versioned model bundle that declares its own operating range.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("predict", help="score one isolate against one model bundle")
    pr.add_argument("--input", required=True, help="ideal-input JSON (GARC variants + covariates)")
    pr.add_argument("--model", required=True, help="model bundle directory")
    pr.add_argument("--threshold", type=float, default=0.5)
    pr.add_argument("--json", action="store_true", help="emit the structured result")
    pr.add_argument("--out", help="write to a file instead of stdout")
    pr.set_defaults(func=_predict)

    de = sub.add_parser("describe", help="print a bundle's declared operating range")
    de.add_argument("--model", required=True, help="model bundle directory")
    de.set_defaults(func=_describe)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
