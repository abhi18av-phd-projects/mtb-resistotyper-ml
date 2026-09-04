"""Layer 1: the curated WHO catalogue, consulted before any model.

Until now the service asserted its own place in the stack — every result carried
``catalogue_consulted: false`` and a caveat saying the models were *intended* for
variants a catalogue grades Unknown or Fail. Intent is not architecture. This
module makes the claim true: the catalogue is consulted first, and a model is
invoked only for the variants the catalogue could not grade.

The catalogue is WHO v2.1 in GARC1, from oxfordmmm/tuberculosis_amr_catalogues,
read by piezo. Two properties earn it that role over a table maintained here.
It is the reference implementation the GARC nomenclature is defined against, so
a mutation name that means one thing to the catalogue means the same thing to
the feature schema. And it is versioned by file name, so which edition graded a
variant is recorded rather than remembered.

The layers do not vote. A catalogue call is returned as the answer and the model
is not consulted for that drug at all, because a curated assessment backed by
phenotype evidence is not improved by averaging it with a logistic regression.
The model exists for the variants where the catalogue is silent, and nowhere
else.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

CATALOGUE = Path(os.environ.get(
    "MTB_CATALOGUE", "/opt/catalogue/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv"))

# piezo's prediction letters. F is the catalogue's explicit "this variant fails
# to be gradable", which is a different statement from U, "this variant is not
# in the catalogue at all". Both are the models' remit; neither is a call.
GRADED = {"R": "R", "S": "S"}
UNGRADED = {"U", "F"}


class CatalogueUnavailable(RuntimeError):
    """piezo or the catalogue file is missing, so Layer 1 cannot run."""


_RESISTANCE = None


def _load():
    global _RESISTANCE
    if _RESISTANCE is not None:
        return _RESISTANCE
    try:
        import piezo
    except ImportError as exc:                      # pragma: no cover
        raise CatalogueUnavailable("piezo is not installed") from exc
    if not CATALOGUE.exists():
        raise CatalogueUnavailable(f"no catalogue at {CATALOGUE}")
    _RESISTANCE = piezo.ResistanceCatalogue(str(CATALOGUE))
    return _RESISTANCE


def version() -> dict:
    """What the catalogue says about itself, for the provenance record.

    Read from the CSV's own columns rather than parsed out of the file name. The
    Oxford naming convention is informative but not load-bearing: a catalogue
    downloaded to ``who21.csv`` is the same catalogue, and a provenance record
    that goes blank because someone renamed a file is worse than useless.
    """
    if not CATALOGUE.exists():
        return {"file": CATALOGUE.name, "error": "not found"}
    with CATALOGUE.open() as fh:
        row = next(csv.DictReader(fh), {}) or {}
    return {"file": CATALOGUE.name,
            "catalogue": row.get("CATALOGUE_NAME"),
            "version": row.get("CATALOGUE_VERSION"),
            "grammar": row.get("CATALOGUE_GRAMMAR"),
            "prediction_values": row.get("PREDICTION_VALUES"),
            "reference": row.get("GENBANK_REFERENCE")}


def grade(variants: list, drug: str) -> dict:
    """Grade one isolate's variants against the catalogue for one drug.

    Returns the catalogue's call and the mutations that produced it. A drug with
    no graded mutation is returned as ``None``, which is the signal to fall
    through to the model.
    """
    cat = _load()
    hits, ungraded = [], []
    for v in variants:
        # The ideal-input instance carries plain dicts; vcf_to_garc yields Variant
        # dataclasses. Both name the same thing, so accept either rather than
        # forcing the caller to know which side of the pipeline it is on.
        gene = v["gene"] if isinstance(v, dict) else v.gene
        mut = v["mutation"] if isinstance(v, dict) else v.mutation
        garc = f"{gene}@{mut}"
        try:
            pred = cat.predict(garc, verbose=False)
        except Exception:
            # A mutation piezo cannot parse is ungraded, not susceptible. Saying
            # S here would turn a parse failure into a clinical-shaped negative.
            ungraded.append(garc)
            continue
        # piezo answers in two shapes and they mean different things. A dict
        # keyed by drug is a per-drug assessment drawn from a catalogue row. A
        # bare string -- almost always "S" -- is the catalogue's DEFAULT RULE
        # firing because nothing matched: it says this mutation is not a known
        # determinant, which is not the same as a curated judgement that this
        # isolate is susceptible to this drug.
        #
        # Treating the bare "S" as a Layer 1 grade is what the layer tests caught:
        # embA@c-12t and PPE69@T19K both return bare "S", so the catalogue
        # swallowed every drug and the models never fired at all. The whole
        # premise inverts. Only a dict is a grade.
        if not isinstance(pred, dict):
            ungraded.append(garc)
            continue
        call = pred.get(drug)
        if call in GRADED:
            hits.append({"mutation": garc, "call": GRADED[call]})
        elif call in UNGRADED or call is None:
            ungraded.append(garc)

    if any(h["call"] == "R" for h in hits):
        # Any single graded resistance mutation confers resistance. This is the
        # catalogue's own logic, not an aggregation choice made here.
        return {"prediction": "R", "layer": 1, "graded": hits, "ungraded": ungraded,
                "catalogue": version()}
    if hits:
        return {"prediction": "S", "layer": 1, "graded": hits, "ungraded": ungraded,
                "catalogue": version()}
    return {"prediction": None, "layer": 1, "graded": [], "ungraded": ungraded,
            "catalogue": version()}
