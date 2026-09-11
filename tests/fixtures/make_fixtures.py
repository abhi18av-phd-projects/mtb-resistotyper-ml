#!/usr/bin/env python3
"""Generate the synthetic VCF fixtures used by CI.

Written rather than captured. A fixture cut from a real isolate carries that
isolate's genotype, and a public test corpus is the wrong place for one. These
files contain nothing but coordinates and alleles already published in the WHO
catalogue, under invented sample names.

They deliberately reproduce the three shapes that have broken this tool before,
because a fixture that exercises only the easy path is a fixture that passes
while the tool is broken:

  * a contig named ``NC-000962-3-H37Rv`` rather than the accession, which is what
    MAGMA emits and what the normaliser has to alias;
  * ``<NON_REF>`` symbolic alleles from gVCF reference-confidence blocks, which
    must be dropped only AFTER multiallelics are split;
  * haploid ``GT`` (``1``, not ``1/1``), which is correct for a haploid organism
    and which gumpy 1.3.8 crashes on unless ploidy is fixed first.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
CONTIG_MAGMA = "NC-000962-3-H37Rv"
CONTIG_PLAIN = "NC_000962.3"
LENGTH = 4411532

# Coordinates and alleles only; both are catalogue-level public facts.
#
# Every REF here is the base H37Rv actually carries at that position, in GENOME
# orientation. That is not a detail: three of these were wrong, and the failure
# is silent. gumpy compares REF against the reference, finds a mismatch, and
# derives no mutation -- so a fixture named for katG S315T produced "INH
# susceptible, 0 mutations carried" and still passed every parser test, because
# the parser tests never look at what the variant means. Run this generator with
# --verify to have it check itself against the reference.
#
# katG is on the MINUS strand, so its genome-orientation base is the complement
# of the one a codon table suggests: S315T is the gene-level AGC->ACC, which on
# the genome reads C->G at 2155168, not G->C.
RPOB_RRDR = [(761139, "C", "T"), (761140, "A", "G")]   # rifampicin, rpoB H445C
KATG_315  = [(2155168, "C", "G")]                       # isoniazid, katG S315T (minus strand)
GYRA_94   = [(7582, "A", "C")]                          # fluoroquinolone, gyrA D94A
OFF_PANEL = [(1000000, "T", "C")]                       # inside gltA2: in no model schema
                                                        # and no catalogue rule, so nothing reads it

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={contig},length={length}>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


def _row(contig, pos, ref, alt, calls, *, non_ref=True, haploid=True):
    alts = f"{alt},<NON_REF>" if non_ref else alt
    fields = [contig, str(pos), ".", ref, alts, "500", "PASS", "DP=120",
              "GT:AD:DP:GQ"]
    for carries in calls:
        gt = ("1" if carries else "0") if haploid else ("1/1" if carries else "0/0")
        ad = "2,118,0" if non_ref else "2,118"
        fields.append(f"{gt}:{ad}:120:99")
    return "\t".join(fields)


def write(name, samples, carriers, *, contig=CONTIG_MAGMA, non_ref=True, haploid=True):
    """`carriers` is one bool per sample per RRDR row."""
    body = [HEADER.format(contig=contig, length=LENGTH, samples="\t".join(samples))]
    for (pos, ref, alt), calls in zip(RPOB_RRDR, carriers, strict=True):
        body.append(_row(contig, pos, ref, alt, calls, non_ref=non_ref, haploid=haploid) + "\n")
    (HERE / name).write_text("".join(body))
    return HERE / name


def write_rows(name, samples, rows, *, contig=CONTIG_MAGMA, non_ref=True,
               haploid=True, header_only=False):
    """Explicit (position, ref, alt, [carrier flags]) rows."""
    body = [HEADER.format(contig=contig, length=LENGTH, samples="\t".join(samples))]
    if not header_only:
        for pos, ref, alt, calls in rows:
            body.append(_row(contig, pos, ref, alt, calls,
                             non_ref=non_ref, haploid=haploid) + "\n")
    (HERE / name).write_text("".join(body))
    return HERE / name


def _verify(gbk: Path) -> None:
    """Check every emitted REF against H37Rv, and fail loudly if one is wrong.

    Optional because the reference is a 10 MB download and CI does not need it
    to test parsing. Worth running whenever a coordinate changes: a wrong REF
    costs nothing at generation time and everything at interpretation time.

        python make_fixtures.py --verify NC_000962.3.gbk
    """
    from Bio import SeqIO

    seq = next(SeqIO.parse(str(gbk), "genbank")).seq
    bad = []
    for path in sorted(HERE.glob("synthetic_*.vcf")):
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            c = line.split("\t")
            pos, ref = int(c[1]), c[3]
            actual = str(seq[pos - 1:pos - 1 + len(ref)])
            if ref != actual:
                bad.append(f"{path.name}: {pos} REF={ref} but H37Rv reads {actual}")
    for b in bad:
        print("MISMATCH", b)
    print(f"verified against {gbk.name}: "
          f"{'all REF bases match' if not bad else f'{len(bad)} MISMATCHED'}")
    if bad:
        raise SystemExit(1)


if __name__ == "__main__":
    # 1. single sample, the full MAGMA shape: dashed contig, NON_REF, haploid GT
    write("synthetic_single.g.vcf", ["SYNTH01"], [[True], [True]])
    # 2. a cohort, so the per-sample split is exercised. Only the first carries.
    write("synthetic_cohort3.g.vcf", ["SYNTH01", "SYNTH02", "SYNTH03"],
          [[True, False, False], [True, False, False]])
    # 3. the plain shape: accession contig, no NON_REF, diploid GT. If this one
    #    fails and the others pass, the fault is in the conversion, not the
    #    normaliser's aliasing and ploidy handling.
    write("synthetic_plain.vcf", ["SYNTH01"], [[True], [True]],
          contig=CONTIG_PLAIN, non_ref=False, haploid=False)
    def one(v):
        return [v]

    # --- drug and layer coverage -------------------------------------------
    write_rows("synthetic_inh_katg.g.vcf", ["SYNTH04"],
               [(p, r, a, one(True)) for p, r, a in KATG_315])
    write_rows("synthetic_mdr.g.vcf", ["SYNTH05"],
               [(p, r, a, one(True)) for p, r, a in RPOB_RRDR + KATG_315])
    write_rows("synthetic_fq_gyra.g.vcf", ["SYNTH06"],
               [(p, r, a, one(True)) for p, r, a in GYRA_94])
    # No resistance mutation anywhere: every drug should fall to Layer 2.
    write_rows("synthetic_susceptible.g.vcf", ["SYNTH07"],
               [(p, r, a, one(False)) for p, r, a in RPOB_RRDR])
    # A variant nothing in the panel reads. It must not change any call.
    write_rows("synthetic_off_panel.g.vcf", ["SYNTH08"],
               [(p, r, a, one(True)) for p, r, a in OFF_PANEL])

    # --- input shapes the normaliser has to survive -------------------------
    # Two real ALTs at one position: multiallelics must be split BEFORE the
    # NON_REF filter, or excluding NON_REF deletes the whole record.
    write_rows("synthetic_multiallelic.g.vcf", ["SYNTH09"],
               [(761139, "C", "T,A", one(True))])
    write_rows("synthetic_insertion.g.vcf", ["SYNTH10"],
               [(761139, "C", "CTT", one(True))])
    # H37Rv reads CAC at 761139-761141, so a two-base deletion anchored there
    # is CAC->C. "CTT" was not the reference and the record was unparseable.
    write_rows("synthetic_deletion.g.vcf", ["SYNTH11"],
               [(761139, "CAC", "C", one(True))])
    # A reference-confidence block and nothing else.
    write_rows("synthetic_nonref_only.g.vcf", ["SYNTH12"],
               [(761139, "C", "<NON_REF>", one(False))], non_ref=False)
    # Diploid GT on a haploid organism: wrong, but callers emit it.
    write_rows("synthetic_diploid_gt.g.vcf", ["SYNTH13"],
               [(p, r, a, one(True)) for p, r, a in RPOB_RRDR], haploid=False)
    # A third contig spelling, neither the accession nor MAGMA's.
    write_rows("synthetic_contig_alias.g.vcf", ["SYNTH14"],
               [(p, r, a, one(True)) for p, r, a in RPOB_RRDR], contig="Chromosome")
    # Header, no records. Must fail with an explanation, not a stack trace.
    write_rows("synthetic_empty.g.vcf", ["SYNTH15"], [], header_only=True)
    # A large-ish cohort, to exercise the per-sample split beyond three.
    write_rows("synthetic_cohort8.g.vcf", [f"SYNTH{i:02d}" for i in range(20, 28)],
               [(p, r, a, [i % 2 == 0 for i in range(8)]) for p, r, a in RPOB_RRDR])

    for f in sorted(HERE.glob("synthetic_*.vcf")):
        print(f"{f.name:36s} {f.stat().st_size:>6d} bytes")

    if "--verify" in sys.argv:
        gbk = sys.argv[sys.argv.index("--verify") + 1]
        _verify(Path(gbk))
