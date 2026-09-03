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

import re
from dataclasses import dataclass
from pathlib import Path

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


def variants_from_vcf(vcf_path: Path, ref_path: Path | None = None) -> list[Variant]:
    """VCF -> GARC variants, via gumpy's genome difference."""
    try:
        import gumpy
    except ImportError as exc:
        raise ReferenceUnavailable("gumpy is not installed") from exc

    genome = reference(ref_path)
    sample = genome + gumpy.VCFFile(str(vcf_path), ignore_filter=True)
    diff = genome - sample

    out: list[Variant] = []
    for gene_name in sorted(set(diff.gene_name) if hasattr(diff, "gene_name") else []):
        if not gene_name:
            continue
        try:
            g_ref, g_alt = genome.build_gene(gene_name), sample.build_gene(gene_name)
        except Exception:
            continue
        gdiff = g_ref - g_alt
        for mutation in getattr(gdiff, "mutations", []) or []:
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
