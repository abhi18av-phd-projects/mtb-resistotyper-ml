"""Generate the demo VCFs the public deployment lists in its object dropdown.

A public prediction service with an empty file list teaches nobody anything, and
the group bucket's real VCFs cannot be the demo: their keys carry sample
identifiers, and enumeration of that prefix would disclose who was sequenced.
So the demo objects are synthetic, and they are generated here rather than
hand-written.

Hand-written was tried and was wrong. The repository's `tests/fixtures/*.g.vcf`
carry the *gene-strand* base as REF for genes on the minus strand -- katG S315T
appears as `2155168 G>C` where H37Rv actually reads C at that position -- so
gumpy sees a reference mismatch and derives no mutation at all. The file looks
like a resistant isolate and predicts susceptible. Those fixtures exercise the
VCF parser, which is what they were written for; they do not exercise the
models. Deriving the coordinates from the reference removes the class of error
rather than fixing one instance of it.

Every emitted variant is checked against H37Rv before it is written: the codon
must translate to the wild-type residue named in the request, and the
substitution must produce the mutant one. A mutation that cannot be reached by a
single-nucleotide change is refused rather than approximated with two.

Usage:
    python make_demo_vcfs.py --reference NC_000962.3.gbk --out .
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

CONTIG = "NC-000962-3-H37Rv"


class Reference:
    def __init__(self, gbk: Path):
        rec = next(SeqIO.parse(str(gbk), "genbank"))
        self.seq = rec.seq
        self.length = len(rec.seq)
        self.genes: dict[str, tuple[int, int, int]] = {}
        for f in rec.features:
            if f.type not in ("CDS", "rRNA"):
                continue
            name = (f.qualifiers.get("gene") or [""])[0]
            if name:
                self.genes[name] = (int(f.location.start), int(f.location.end),
                                    f.location.strand)

    def substitution(self, gene: str, wt: str, codon: int, mut: str) -> tuple[int, str, str]:
        """Genome-level (1-based pos, ref, alt) for an amino-acid substitution.

        Returned in GENOME orientation, which is the only orientation a VCF has.
        For a minus-strand gene the base written is the complement of the one the
        codon table would suggest, and getting that backwards is precisely the
        bug this function exists to prevent.
        """
        start, end, strand = self.genes[gene]
        cds = self.seq[start:end]
        if strand == -1:
            cds = cds.reverse_complement()
        i = (codon - 1) * 3
        ref_codon = cds[i:i + 3]
        observed = str(Seq(str(ref_codon)).translate())
        if observed != wt:
            raise ValueError(
                f"{gene} codon {codon} is {observed}, not {wt} — refusing to emit")

        for offset in range(3):
            for base in "ACGT":
                if base == str(ref_codon[offset]):
                    continue
                alt_codon = ref_codon[:offset] + base + ref_codon[offset + 1:]
                if str(Seq(str(alt_codon)).translate()) != mut:
                    continue
                gene_index = i + offset
                if strand == 1:
                    pos0 = start + gene_index
                    return pos0 + 1, str(self.seq[pos0]), base
                pos0 = end - 1 - gene_index
                return (pos0 + 1, str(self.seq[pos0]),
                        str(Seq(base).complement()))
        raise ValueError(
            f"{gene} {wt}{codon}{mut} needs more than one nucleotide change")

    def promoter(self, gene: str, offset: int, wt: str, mut: str) -> tuple[int, str, str]:
        """Genome-level call for a `c-15t`-style promoter variant.

        Offset is negative and counted in the GENE's direction from the first
        base of the start codon, which is how GARC names it.
        """
        start, end, strand = self.genes[gene]
        if strand == 1:
            pos0 = start + offset
            ref, alt = wt.upper(), mut.upper()
        else:
            pos0 = end - 1 - offset
            ref = str(Seq(wt.upper()).complement())
            alt = str(Seq(mut.upper()).complement())
        actual = str(self.seq[pos0])
        if actual != ref:
            raise ValueError(
                f"{gene} {offset} reads {actual}, not {ref} — refusing to emit")
        return pos0 + 1, ref, alt


HEADER = f"""##fileformat=VCFv4.2
##source=mtb-resistotyper-ml demo generator (synthetic; not a real isolate)
##reference=NC_000962.3
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={CONTIG},length=4411532>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"""


def write_vcf(path: Path, samples: list[str],
              calls: dict[str, list[tuple[int, str, str]]]) -> None:
    """One VCF, one or many samples. `calls` maps sample id -> its variants."""
    rows: dict[tuple[int, str, str], set[str]] = {}
    for sample, variants in calls.items():
        for v in variants:
            rows.setdefault(v, set()).add(sample)

    lines = [HEADER, "#" + "\t".join(
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples)]
    for (pos, ref, alt) in sorted(rows):
        carriers = rows[(pos, ref, alt)]
        gts = ["1:2,118,0:120:99" if s in carriers else "0:118,2,0:120:99"
               for s in samples]
        lines.append("\t".join([CONTIG, str(pos), ".", ref, f"{alt},<NON_REF>",
                                "500", "PASS", "DP=120", "GT:AD:DP:GQ"] + gts))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    a = ap.parse_args()
    ref = Reference(a.reference)
    a.out.mkdir(parents=True, exist_ok=True)

    katG = ref.substitution("katG", "S", 315, "T")
    rpoB = ref.substitution("rpoB", "S", 450, "L")
    gyrA = ref.substitution("gyrA", "D", 94, "A")
    embB = ref.substitution("embB", "M", 306, "V")
    fabG1 = ref.promoter("fabG1", -15, "c", "t")
    print("katG S315T ", katG)
    print("rpoB S450L ", rpoB)
    print("gyrA D94A  ", gyrA)
    print("embB M306V ", embB)
    print("fabG1 c-15t", fabG1)

    # A susceptible isolate needs a variant that is real but consequence-free,
    # otherwise it tests nothing: an empty VCF and a wild-type one travel
    # different code paths, and it is the wild-type one a user will bring.
    katG_syn = ref.substitution("katG", "P", 365, "P")

    write_vcf(a.out / "demo-01-susceptible.vcf", ["DEMO01"], {"DEMO01": [katG_syn]})
    write_vcf(a.out / "demo-02-inh-katG-S315T.vcf", ["DEMO02"], {"DEMO02": [katG]})
    write_vcf(a.out / "demo-03-rif-rpoB-S450L.vcf", ["DEMO03"], {"DEMO03": [rpoB]})
    write_vcf(a.out / "demo-04-mdr-katG-rpoB.vcf", ["DEMO04"], {"DEMO04": [katG, rpoB]})
    write_vcf(a.out / "demo-05-pre-xdr-mdr-gyrA-D94A.vcf", ["DEMO05"],
              {"DEMO05": [katG, rpoB, gyrA]})
    write_vcf(a.out / "demo-06-inh-promoter-fabG1-c-15t.vcf", ["DEMO06"],
              {"DEMO06": [fabG1]})
    write_vcf(a.out / "demo-07-cohort-of-four.vcf",
              ["DEMO07", "DEMO08", "DEMO09", "DEMO10"],
              {"DEMO07": [katG_syn], "DEMO08": [katG], "DEMO09": [rpoB, embB],
               "DEMO10": [katG, rpoB, gyrA]})
    for f in sorted(a.out.glob("demo-*.vcf")):
        print("wrote", f.name, f.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
