"""Turn a VCF into the GARC variants the deployed models are keyed on.

The models' features are amino-acid and promoter changes named in GARC, for
example `rpoB_S450L`. Deriving those from a VCF is codon-level annotation
against H37Rv, which is exactly what gnomonicus does, so this delegates rather
than reimplementing it badly.

Two consequences worth stating rather than discovering. A variant with no gene
cannot be expressed in GARC at all, so intergenic calls are dropped here and are
invisible to the model however important they may be. And when the reference is
unavailable the conversion is refused rather than approximated: a wrong mutation
name silently scores as a reference allele, which looks like a confident
susceptible call.
"""

from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

from mtb_resistotyper_ml import catalogue as _catalogue

REFERENCE_ENV = "MTB_REFERENCE_GENBANK"


@dataclass(frozen=True)
class Variant:
    gene: str
    mutation: str


class ReferenceUnavailable(RuntimeError):
    """Raised when GARC cannot be derived, rather than guessed at."""


def _load_reference(path: Path):
    try:
        import gumpy
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ReferenceUnavailable(
            "gumpy is not installed, so a VCF cannot be converted to GARC. "
            "Submit variants in the ideal-input JSON format instead."
        ) from exc
    if not path.exists():
        raise ReferenceUnavailable(
            f"no H37Rv reference at {path}. Set {REFERENCE_ENV} to a GenBank file."
        )
    return gumpy.Genome(str(path))


_CACHE: dict[str, object] = {}


def reference(path: Path | None = None):
    p = Path(path or Path(__file__).with_name("NC_000962.3.gbk"))
    key = str(p)
    if key not in _CACHE:
        _CACHE[key] = _load_reference(p)
    return _CACHE[key]


def genes_of_interest(models_dir: Path | None = None) -> set[str]:
    """Genes any layer can actually act on.

    A clinical isolate differs from H37Rv in hundreds of genes, and gumpy's
    per-gene diff costs roughly a fifth of a second each. Building all of them
    took two to four minutes per sample, nearly all of it spent deriving PE and
    PE_PGRS mutations that no model has a feature for and no catalogue entry
    mentions. Restricting to the union of the model schemas and the catalogue's
    own gene set is not an optimisation that trades away information: a gene
    outside that union cannot change any answer this service gives.
    """
    genes: set[str] = set()
    models = Path(models_dir or os.environ.get("MTB_MODELS", "/opt/models"))
    if models.exists():
        # Both layouts, because a deployment carrying two compendium releases
        # nests them one level deeper (MODELS/<release>/<DRUG>/) and a gene
        # missed here is not an error anywhere -- it is silently never
        # reconstructed, so the variant that needed it scores as wild type.
        schemas = list(models.glob("*/feature_schema.json"))
        schemas += list(models.glob("*/*/feature_schema.json"))
        for schema in schemas:
            for f in json.loads(schema.read_text()).get("features", []):
                if f.get("gene"):
                    genes.add(f["gene"])
    # Reads catalogue.CATALOGUE rather than re-deriving MTB_CATALOGUE here: this
    # is the same path --catalogue sets (the CLI writes it into that module, not
    # into the environment) and the same default catalogue.py itself falls back
    # to. Re-deriving it independently, from `os.environ.get("MTB_CATALOGUE",
    # "")`, was a second bug on top of missing --catalogue: an unset env var
    # made this `Path("")`, which is `Path(".")`, which exists, so the code
    # went on to open the current directory as a file.
    cat = _catalogue.CATALOGUE
    if cat.exists():
        with cat.open() as fh:
            for row in csv.DictReader(fh):
                m = row.get("MUTATION", "")
                if "@" in m:
                    genes.add(m.split("@", 1)[0])
    return genes


def variants_from_vcf(vcf_path: Path, ref_path: Path | None = None,
                      genes: set[str] | None = None) -> list[Variant]:
    """VCF -> GARC variants, via gumpy's genome difference.

    ``genes`` restricts which genes are reconstructed; it defaults to every gene
    the models or the catalogue can act on. Pass an empty set to disable the
    restriction and reconstruct everything, which is slow and only useful when
    exploring what a sample carries outside the current feature space.
    """
    try:
        import gumpy
    except ImportError as exc:
        raise ReferenceUnavailable("gumpy is not installed") from exc

    genome = reference(ref_path)
    sample = genome + gumpy.VCFFile(str(vcf_path), ignore_filter=True)
    diff = genome - sample

    out: list[Variant] = []
    # gene_name carries None for intergenic differences, and numpy string arrays
    # do not sort against None. Intergenic positions are dropped rather than
    # coerced: this feature space is gene-anchored GARC, so a variant with no
    # gene has no name in it.
    names = {str(g) for g in getattr(diff, "gene_name", []) if g is not None and str(g)}
    wanted = genes_of_interest() if genes is None else genes
    if wanted:
        names &= wanted
    for gene_name in sorted(names):
        try:
            g_ref, g_alt = genome.build_gene(gene_name), sample.build_gene(gene_name)
        except Exception:
            continue
        gdiff = g_ref - g_alt
        # `or []` on a numpy array raises "truth value ... is ambiguous", and
        # gumpy returns arrays here. Test for None explicitly rather than
        # leaning on truthiness, which numpy deliberately refuses to define.
        muts = getattr(gdiff, "mutations", None)
        for mutation in ([] if muts is None else list(muts)):
            out.append(Variant(gene=gene_name, mutation=str(mutation)))
    return out


_SAMPLE_LINE = re.compile(r"^#CHROM\s")


def sample_ids(vcf_path: Path) -> list[str]:
    """Sample columns in a VCF header, so a multi-sample file names its members."""
    for line in Path(vcf_path).read_text(errors="replace").splitlines():
        if _SAMPLE_LINE.match(line):
            cols = line.split("\t")
            return cols[9:] if len(cols) > 9 else []
        if not line.startswith("#"):
            break
    return []


def to_instance(sample_id: str, variants: list[Variant],
                lineage: str = "unknown", median_coverage: float | None = None,
                breadth: float | None = None) -> dict:
    """Build the ideal-input instance the runner's feature schema expects."""
    return {
        "schema_version": "1.0.0",
        "sample_id": sample_id,
        "nomenclature": "GARC",
        "source": {"producer": "gumpy/gnomonicus", "reference": "NC_000962.3"},
        "covariates": {
            "lineage": lineage,
            # Absent covariates are mean-imputed by the scorer, which is the
            # documented behaviour: a genotype-only call should not be swamped
            # by an unknown quality covariate.
            **({"median_coverage": median_coverage} if median_coverage is not None else {}),
            **({"breadth": breadth} if breadth is not None else {}),
        },
        "variants": [{"gene": v.gene, "mutation": v.mutation} for v in variants],
    }
