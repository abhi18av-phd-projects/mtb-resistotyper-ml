"""Make a real caller's VCF acceptable to gumpy, using bcftools throughout.

Written against genuine MAGMA output pulled from the cluster, which failed the
converter in four distinct ways that no hand-written VCF reproduces. Each is a
documented bcftools operation; none is parsed here by hand.

1. Contig naming. MAGMA writes ``NC-000962-3-H37Rv``; the GenBank record is
   ``NC_000962.3``. ``bcftools annotate --rename-chrs``. The alias table is the
   one thing bcftools cannot know, since only this project knows those spellings
   mean H37Rv.

2. Multiallelic sites. ``bcftools norm -m -any``. This must come FIRST, because
   it is what makes step 3 correct rather than total.

3. gVCF reference blocks. HaplotypeCaller with ``--emit-ref-confidence GVCF``
   puts ``<NON_REF>`` on every row, so excluding it before splitting deletes the
   entire file — verified: 0 records survived. After splitting, a record whose
   single ALT is ``<NON_REF>`` is a reference block and only those are dropped.

4. Haploid genotypes. GATK ``--sample-ploidy 1`` is correct for a haploid
   organism, and gumpy 1.3.8 crashes on it: ``difference.py:268`` evaluates
   ``gt[0] == gt[1]`` on a one-element tuple. ``bcftools +fixploidy`` rewrites
   ``1`` as ``1/1``. This is a workaround for an upstream defect, not a
   correction to the data, and it should be removed when gumpy handles haploid
   calls.

Multi-sample cohorts are then split one file per sample with ``bcftools +split``,
dropping hom-ref and no-call: a site where a sample has no variant carries no
information about it, and emitting it would assert wild type at a position that
may simply not have been called.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

ALIASES = {
    "NC-000962-3-H37Rv": "NC_000962.3",   # MAGMA / XBS
    "NC_000962": "NC_000962.3",
    "NC_000962_3": "NC_000962.3",
    "AL123456.3": "NC_000962.3",          # the EMBL accession for the same sequence
    "AL123456": "NC_000962.3",
    "H37Rv": "NC_000962.3",
    "Chromosome": "NC_000962.3",
}
# MTBseq calls against a reconstructed ancestral genome that differs from H37Rv
# at roughly a thousand positions. Renaming it would shift every coordinate, so
# it is refused rather than aliased.
REFUSED = {"MTB_anc": "a reconstructed ancestral genome, not H37Rv"}


class NormalisationFailed(RuntimeError):
    pass


def available() -> bool:
    return shutil.which("bcftools") is not None


def _run(args: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(args, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise NormalisationFailed(f"{' '.join(args[:3])}…: {r.stderr.strip()[:400]}")
    return r


def contigs(vcf: Path) -> list[str]:
    r = _run(["bcftools", "query", "-f", "%CHROM\n", str(vcf)])
    seen, out = set(), []
    for c in r.stdout.split():
        if c not in seen:
            seen.add(c); out.append(c)
    return out


def normalise(vcf: Path, workdir: Path, *, genes: set[str] | None = None,
              genbank: Path | None = None) -> list[tuple[str, Path]]:
    """Return ``[(sample_name, single_sample_vcf)]`` ready for gumpy.

    When ``genes`` and ``genbank`` are supplied the records are cut down to those
    genes first, which is what makes the conversion fast enough to serve
    synchronously.
    """
    if not available():
        raise NormalisationFailed("bcftools is not installed")
    workdir.mkdir(parents=True, exist_ok=True)

    for c in contigs(vcf):
        if c in REFUSED:
            raise NormalisationFailed(
                f"contig {c!r} is {REFUSED[c]}; its coordinates differ from H37Rv, so "
                f"mapping it would relabel real mutations. Re-call against H37Rv.")
        if c not in ALIASES and c != "NC_000962.3":
            raise NormalisationFailed(
                f"contig {c!r} is not a recognised spelling of H37Rv (NC_000962.3). "
                f"Refusing rather than assuming the reference.")

    chrmap = workdir / "chrmap.txt"
    chrmap.write_text("".join(f"{k}\t{v}\n" for k, v in ALIASES.items()))
    ploidy = workdir / "ploidy.txt"
    ploidy.write_text("* * * * 2\n")

    clean = workdir / "clean.vcf.gz"
    p1 = subprocess.Popen(["bcftools", "annotate", "--rename-chrs", str(chrmap), str(vcf), "-Ou"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["bcftools", "norm", "-m", "-any", "-Ou"],
                          stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p3 = subprocess.Popen(["bcftools", "view", "-e", 'ALT="<NON_REF>" || ALT="."', "-Ou"],
                          stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p4 = subprocess.Popen(["bcftools", "view", "--trim-alt-alleles", "-Oz", "-o", str(clean)],
                          stdin=p3.stdout, stderr=subprocess.PIPE)
    for p in (p1, p2, p3):
        p.stdout.close()
    err = p4.communicate()[1]
    if p4.returncode != 0:
        raise NormalisationFailed(f"bcftools clean-up: {err.strip()[:400]}")

    if genes and genbank:
        bed = regions_bed(genes, genbank, workdir / "regions.bed")
        if bed is not None:
            restricted = workdir / "restricted.vcf.gz"
            clean = restrict(clean, bed, restricted)

    samples = _run(["bcftools", "query", "-l", str(clean)]).stdout.split()
    if not samples:
        raise NormalisationFailed("the VCF carries no sample column, so nothing can be genotyped")

    _run(["bcftools", "index", "-f", str(clean)])
    split = workdir / "split"
    split.mkdir(exist_ok=True)
    _run(["bcftools", "+split", str(clean), "-i", 'GT!="ref" && GT!="mis"',
          "-Ov", "-o", str(split)])

    out = []
    for s in samples:
        raw = split / f"{s}.vcf"
        if not raw.exists():
            continue
        fixed = workdir / f"{s}.ready.vcf"
        sfile = workdir / f"{s}.sex.txt"
        sfile.write_text(f"{s}\tM\n")
        r = subprocess.run(["bcftools", "+fixploidy", str(raw), "--",
                            "-p", str(ploidy), "-s", str(sfile)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise NormalisationFailed(f"bcftools +fixploidy: {r.stderr.strip()[:300]}")
        fixed.write_text(r.stdout)
        out.append((s, fixed))
    if not out:
        raise NormalisationFailed(
            "no sample carried a non-reference call after normalisation")
    return out


def regions_bed(genes: set[str], genbank: Path, dest: Path) -> Path | None:
    """A BED of the genes any layer can act on, for restricting the VCF.

    This is where the real cost sits. gumpy applies a VCF to the genome at
    roughly five variants a second, so a clinical isolate's twelve hundred calls
    take four minutes before a single gene is reconstructed. Restricting the gene
    loop afterwards does not help: the work has already been done.

    Sixty genes carry every feature the models use and every gene the catalogue
    names. A variant outside them cannot change any answer this service gives, so
    cutting them out of the VCF first is not a speed-versus-accuracy trade; it is
    declining to compute a result nothing reads. Coordinates come from the same
    GenBank record gumpy parses, so the two cannot disagree about where a gene is.
    """
    import re as _re
    if not genes or not genbank.exists():
        return None
    # Parse the feature table directly. Biopython would do this too, but the
    # image already reads this file once through gumpy and a second full parse
    # costs more than the scan it replaces.
    rows: list[tuple[int, int, str]] = []
    loc = None
    text = genbank.read_text(errors="replace")
    for block in text.split("\n     gene            ")[1:]:
        head, _, rest = block.partition("\n")
        m = _re.search(r"(\d+)\.\.(\d+)", head)
        if not m:
            continue
        name = None
        for line in rest.splitlines():
            g = _re.search(r'/(?:gene|locus_tag)="([^"]+)"', line)
            if g:
                name = g.group(1)
                break
            if line.startswith("     ") and not line.startswith("                     "):
                break
        if name and name in genes:
            # BED is half-open and zero-based; GenBank is inclusive and one-based.
            # Pad by 100 bp so promoter variants, which GARC names with negative
            # offsets such as fabG1@c-15t, are not cut away with the intergenic.
            rows.append((max(0, int(m.group(1)) - 101), int(m.group(2)) + 100, name))
    if not rows:
        return None
    rows.sort()
    dest.write_text("".join(f"NC_000962.3\t{s}\t{e}\t{n}\n" for s, e, n in rows))
    return dest


def restrict(vcf: Path, bed: Path, out: Path) -> Path:
    """Keep only records inside the regions of interest."""
    _run(["bcftools", "index", "-f", str(vcf)])
    _run(["bcftools", "view", "-R", str(bed), str(vcf), "-Oz", "-o", str(out)])
    return out
