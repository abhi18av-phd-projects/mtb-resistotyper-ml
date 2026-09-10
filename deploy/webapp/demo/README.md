# Demo isolates for the public deployment

Seven synthetic VCFs, listed in the service's object dropdown so that a reader
arriving from the paper can run a prediction without having a VCF of their own.
They live in the group bucket under `demo/`, which is the only real-data-free
prefix `MTB_LIST_PREFIXES` enumerates.

**No file here is a real isolate.** That is the point: object keys in a research
group's bucket carry sample identifiers, so the prefix a public deployment is
allowed to enumerate cannot contain sequencing. These are constructed.

## Why they are generated rather than written

`tests/fixtures/*.g.vcf` cannot serve as the demo. For genes on the minus strand
those fixtures carry the **gene-strand base** as REF: katG S315T appears as
`2155168 G>C`, where H37Rv reads C at that position and the true call is `C>G`.
gumpy sees a reference mismatch, derives no mutation, and the file that is named
for a resistance marker predicts **susceptible with zero mutations carried**.
They exercise the VCF parser, which is what they were written for.

`make_demo_vcfs.py` derives every coordinate from the reference instead, and
refuses to emit a variant whose codon does not translate to the residue named.

## Regenerating

    curl -sSfL https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues/master/catalogues/NC_000962.3/NC_000962.3.gbk -o ref.gbk
    python make_demo_vcfs.py --reference ref.gbk --out .
    abc data push "./*.vcf" s3://su-mbhg-bioinformatics/demo/

The printed coordinates are the check worth reading: katG S315T is
`2155168 C>G`, rpoB S450L `761155 C>T`, gyrA D94A `7582 A>C`.

## What each one demonstrates

| file | carries | expected |
|---|---|---|
| `demo-01-susceptible` | katG P365P (synonymous) | all susceptible |
| `demo-02-inh-katG-S315T` | katG S315T | INH resistant |
| `demo-03-rif-rpoB-S450L` | rpoB S450L | RIF resistant |
| `demo-04-mdr-katG-rpoB` | both of the above | MDR |
| `demo-05-pre-xdr-mdr-gyrA-D94A` | + gyrA D94A | MDR + fluoroquinolone |
| `demo-06-inh-promoter-fabG1-c-15t` | fabG1 c-15t | INH + ETH resistant |
| `demo-07-cohort-of-four` | four samples, four profiles | one multi-sample VCF |

The resistant calls are all **Layer 1**: these are canonical markers, so the WHO
catalogue grades them and no model is consulted for that drug. Layer 2 answers
every other drug on the same isolate, which is the division the service exists
to make visible.
