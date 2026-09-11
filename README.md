# mtb-resistotyper-ml

Predict *Mycobacterium tuberculosis* drug resistance from variant calls, using a **versioned
model bundle that declares its own operating range**.

```bash
pip install mtb-resistotyper-ml
```

That bare install scores an isolate you have already converted, and pulls in nothing:
the scoring path is a closed form over the standard library. Two capabilities need
more, and each is an extra so you only carry what you use.

| you want to | install | also needs |
|---|---|---|
| score an already-converted isolate | `mtb-resistotyper-ml` | — |
| Layer 1, the curated WHO catalogue | `mtb-resistotyper-ml[catalogue]` | a catalogue CSV |
| read a VCF / gVCF directly | `mtb-resistotyper-ml[vcf]` | **bcftools on `PATH`**, an H37Rv GenBank |
| everything | `mtb-resistotyper-ml[all]` | both of the above |

**bcftools is a system dependency and pip cannot supply it.** Install it with your
package manager (`apt install bcftools`, `brew install bcftools`, `conda install -c
bioconda bcftools`). Without it the library raises `BcftoolsUnavailable` rather than
producing a wrong answer from an unnormalised file.

Both capabilities are imported lazily, so a bare install imports and runs; ask it to do
something it lacks the dependency for and it says which one is missing.

## A reproducible install: dev and prod (pixi)

`pip install` above still works and needs nothing else. [pixi](https://pixi.sh) is the
route to a full, reproducible environment in one command, because it is the one tool
here that can actually supply **bcftools** — a conda package, not a PyPI one, and the
thing the plain-pip install above has to ask you to find yourself.

```bash
pixi install            # dev: editable install, pytest, ruff, bcftools, gumpy, piezo
pixi install -e prod    # prod: a real (non-editable) install of this checkout, same
                         # deps, no dev tooling — what deploy/cli/Dockerfile builds
```

```bash
pixi run test            # pytest -q
pixi run lint             # ruff check src tests
pixi run predict -- --vcf isolate.vcf --models ./models --tsv
```

`pixi.toml` pins `biopython==1.83`: gumpy reads `SeqFeature.strand`, which Biopython
deprecated in 1.81 and removed in 1.85, and an unpinned resolve reliably picks a newer
one — confirmed by running a real VCF through the unpinned environment, which failed
gene reconstruction on every gene, not just some (`app/Dockerfile` carries the identical
pin, for the identical reason).

**A source edit needs `pixi reinstall` under `prod`, not `pixi install` again.** `prod`
installs this checkout non-editable, and pixi does not notice a source change on a plain
`pixi install` for a local path dependency — `dev`'s editable install always picks edits
up live; `prod` needs telling:

```bash
pixi reinstall -e prod --frozen
```

## Two ways in, one computation

The command line and any Python caller enter through the same function, so a call made
at a terminal and a call made inside a service are the same computation on the same
input.

```bash
# a MAGMA gVCF, single- or multi-sample; a cohort file yields one result set per sample
mtb-resistotyper-ml predict --vcf cohort.g.vcf.gz \
    --reference NC_000962.3.gbk --models ./models --tsv
```

```python
from pathlib import Path
from mtb_resistotyper_ml import instances_from_vcf, resolve_all, catalogue_ready

for instance in instances_from_vcf("cohort.g.vcf.gz", "NC_000962.3.gbk",
                                   models="./models"):
    for row in resolve_all(instance, Path("./models"),
                           use_catalogue=catalogue_ready()):
        print(row["sample_id"], row["drug"], row["prediction"], row["source"])
```

`iter_instances_from_vcf` is the streaming form, for a cohort too large to hold at once.
Every row carries `source` and `layer`, so a caller can always tell whether the
catalogue or a model answered.

## The models are not in here

This is the **runner**. The trained models live in
[`mtb-resistotyper-ml-models`](https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-models),
with their own release cadence, and the two are joined only by the artefact contract:
`feature_schema.json` fixes the input, `model_card.json` carries the operating range.

That separation is the point. A revised WHO catalogue edition, a larger training cohort, or a
model retrained for a different population produces a new version *there* and leaves this tool
untouched — and a third party can retrain and publish a bundle this runner will load without
needing anything from us.

```bash
git clone https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-models
mtb-resistotyper-ml describe --model mtb-resistotyper-ml-models/models/RIF
```

```
drug        RIF
model       standardscaler+l2_logistic
features    41
AUC         0.9455   (nested lineage-leave-one-out cross-validation)
tier        usable
per lineage lineage1=0.912  lineage2=0.964  lineage3=0.948  lineage4=0.958
```

## Predicting

```bash
mtb-resistotyper-ml predict \
  --input isolate.json \
  --model mtb-resistotyper-ml-models/models/RIF
```

```
sample     TEST-RIF-R
drug       RIF
call       R   P(R) = 0.9769
tier       usable  (AUC 0.946, nested lineage-leave-one-out)

reasons    2 mutation(s) carried
  rpoB@S450L           coef +2.437   logit +4.186   on-target  primary driver
  katG@S315T           coef +0.616   logit +0.775   OFF-target co-selected passenger (multi-drug linkage)
```

`--json` emits the structured result for machine consumption.

Note what the second line of reasoning does: `katG S315T` genuinely raises the rifampicin
probability, because in this population it travels with `rpoB` mutations — but it is not a
rifampicin mechanism, and the tool says so rather than presenting it as one.

## Running the container image

`deploy/cli/Dockerfile` packages the tool, bcftools, the H37Rv reference and the WHO
catalogue into one image. It does **not** bake in a model — mount a downloaded
[`mtb-resistotyper-ml-models`](https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-models)
release instead, the same directory `--models` takes above, so the image tag pins the
tool version and the model release pins itself:

```bash
docker run --rm \
  -v "$PWD/models:/models:ro" -v "$PWD/data:/data:ro" -v "$PWD/out:/out" \
  ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml/mtb-resistotyper-ml-cli:vX.Y.Z \
  predict --vcf /data/isolate.vcf --models /models --outdir /out
```

Built and pushed by `.github/workflows/docker-cli.yml` on the same `v*` tag
`release.yml` publishes to PyPI from, so the wheel and the image always name the same
commit. `--models` may also point at one release inside the nested
`MODELS/<release>/<DRUG>/` layout `deploy/webapp/bring-up.sh` produces for more than
one installed release.

## Input

A JSON instance validating against `input_spec.schema.json`: GARC-nomenclature variant calls
plus sample covariates, and **no catalogue interpretation**. Replacing that interpretation is
what the model does; including it would leak the answer. Absence of a variant means wild type at
that locus.

```json
{
  "schema_version": "1.0.0",
  "sample_id": "ERR1234567",
  "nomenclature": "GARC",
  "covariates": {"lineage": "lineage2", "median_coverage": 55.0, "breadth": 0.98},
  "variants": [{"gene": "rpoB", "mutation": "S450L"}]
}
```

## What this tool does not do

**It is a Layer 2 predictor.** Layer 1, the curated WHO catalogue via `piezo`, is
consulted first (`mtb_resistotyper_ml.catalogue`); a model is invoked only for a drug
the catalogue graded Unknown, Fail, or not at all — roughly 40% of catalogue entries,
plus every novel variant — and a model never overrides a catalogue call. Every result
carries `layer` and `source` so a caller can always tell which one answered. Passing
`--no-catalogue` (or omitting `--catalogue`/`MTB_CATALOGUE` with the `catalogue` extra
not installed) turns Layer 1 off entirely, and every call then becomes a model call.

**It is not clinically validated.** No prospective validation has been done. Research use only.

**It has no model for bedaquiline, clofazimine, linezolid or delamanid.** Those perform at
chance under honest evaluation, so no bundle is published for them.

## Runtime

Scoring is a closed form over the standard library:

    p(R) = sigmoid( ((x - scaler_mean) / scaler_scale) . coef + intercept )

so this package declares no runtime dependencies, and the scoring path can be vendored or
reimplemented in another language. **Training is a different matter** and is not part of this
package: reproducing a model requires the training pipeline in the paper-assets repository,
which needs a JVM (H2O, for the model-class comparison) and a GPU (protein-language-model
features and the deep-network sweep).

## Why the coefficients are the explanation

The deployed model class was not chosen for interpretability. It was chosen because it *won*:
on identical features and folds, a regularised linear model beat gradient-boosted trees and a
tuned deep network under lineage-leave-one-out evaluation, because the extra capacity was
fitting population structure that does not transfer. The model being small and exactly
explainable is a consequence of honest evaluation, not a trade against accuracy.

So the per-feature logit contributions are not an approximation of the model. They are the
model, and they sum to the logit — which `tests/test_serving.py` asserts.

## Tests

```bash
MTB_MODELS=/path/to/mtb-resistotyper-ml-models pytest
```

Fixtures live in `tests/data` with recorded expected outputs in `tests/expected`, so a change in
scoring cannot pass silently.

## Licence

EPL-2.0.
