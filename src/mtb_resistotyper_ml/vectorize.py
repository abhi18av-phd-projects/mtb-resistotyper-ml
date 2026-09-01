"""Ideal input (the intermediate representation) -> the exact model feature row.

This is the bridge between the ingest contract (``input_spec.schema.json``) and
the model contract (``feature_schema.json``). Getting the row byte-identical to
the training mart row is the entire point of the stage; ``roundtrip`` is the gate
that proves it.

Construction rules mirror the mart builder exactly:

* ``raw__<gene>_<mutation>`` is 1 iff the input carries a variant whose verbatim
  ``gene``/``mutation`` concatenation equals the feature name. The concatenation
  is verbatim and unsanitised, so a wildcard call such as ``rpoB S450X`` simply
  does not match ``raw__rpoB_S450L`` and becomes 0. That is the documented
  default: a failed call is treated as wild type.
* ``cov__lineage_<k>`` is 1 iff the input's lineage equals the k it tests.
* coverage covariates take the input value, or the schema's training-median fill
  when absent.
"""

from __future__ import annotations

from typing import Any

_INT_KINDS = {"mutation", "lineage"}


def vectorize(instance: dict, schema: dict) -> dict[str, Any]:
    """Return an ordered ``{feature_name: value}`` row for ``schema``."""
    present = {f"raw__{v['gene']}_{v['mutation']}"
               for v in instance.get("variants", [])}
    covariates = instance.get("covariates", {})
    lineage = covariates.get("lineage")

    row: dict[str, Any] = {}
    for feature in schema["features"]:
        kind, name = feature["kind"], feature["name"]
        if kind == "mutation":
            value: Any = 1 if name in present else 0
        elif kind == "lineage":
            value = 1 if lineage == feature["lineage"] else 0
        elif kind == "coverage":
            raw = covariates.get(feature["source"])
            value = feature["fill"] if raw is None else float(raw)
        else:
            raise ValueError(f"unknown feature kind {kind!r} for {name!r}")
        row[name] = int(value) if kind in _INT_KINDS else float(value)
    return row


def fired(row: dict[str, Any]) -> list[str]:
    """Mutation features carried by this isolate."""
    return [k for k, v in row.items() if k.startswith("raw__") and float(v) != 0.0]
