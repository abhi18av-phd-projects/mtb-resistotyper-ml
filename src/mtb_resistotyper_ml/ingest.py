"""Turning a VCF into scorable instances — the one implementation.

The command line, the web service and any pipeline that imports this package all
enter here. A prediction is only comparable across those surfaces if the genotype
behind it was derived the same way, and the surest way to keep three callers
identical is to give them one function rather than three copies of it.
"""
from __future__ import annotations

import tempfile
from collections.abc import Iterator
from pathlib import Path

from . import normalise
from .vcf_to_garc import ReferenceUnavailable, genes_of_interest, to_instance, variants_from_vcf

__all__ = [
    "BcftoolsUnavailable",
    "ReferenceUnavailable",
    "instances_from_vcf",
    "iter_instances_from_vcf",
]


class BcftoolsUnavailable(RuntimeError):
    """bcftools is not on PATH; a VCF cannot be normalised without it."""


def iter_instances_from_vcf(vcf: str | Path, reference: str | Path, *,
                            models: str | Path | None = None,
                            lineage: str = "unknown",
                            workdir: str | Path | None = None,
                            ) -> Iterator[dict]:
    """Yield one ideal-input instance per sample in ``vcf``.

    A multi-sample cohort yields several. Scoring only the first column and
    naming the rest after the file reports one isolate's genotype under every
    sample's name, so bcftools splits the file and each sample is carried
    separately.

    ``models`` narrows the conversion to the genes the deployed bundles and the
    catalogue actually read. Passing None converts every gene, which is correct
    but slower.
    """
    ref = Path(reference)
    if not ref.exists():
        raise ReferenceUnavailable(f"reference genbank not found: {ref}")
    if not normalise.available():
        raise BcftoolsUnavailable(
            "bcftools is not on PATH; it is required to normalise a VCF")

    genes = genes_of_interest(Path(models) if models else None)

    def _run(work: Path) -> Iterator[dict]:
        for sample, ready in normalise.normalise(Path(vcf), work,
                                                 genes=genes, genbank=ref):
            yield to_instance(sample, variants_from_vcf(ready, ref, genes=genes),
                              lineage=lineage)

    if workdir is not None:
        yield from _run(Path(workdir))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            yield from _run(Path(tmp))


def instances_from_vcf(vcf: str | Path, reference: str | Path, **kw) -> list[dict]:
    """``iter_instances_from_vcf`` as a list, for callers that want the count."""
    return list(iter_instances_from_vcf(vcf, reference, **kw))
