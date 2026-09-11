# The CLI image. Bring your own VCF and your own downloaded model release;
# this image supplies the tool, bcftools, and the H37Rv reference + WHO
# catalogue Layer 1 needs. It does NOT bake in a model -- that is what a
# release from https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-models
# (or a Zenodo deposit built by package_bundles.py) is for, and baking one in
# would tie an image tag to a model edition instead of to the tool version.
#
# Usage, once you have a model release's `models/` directory locally:
#
#   docker run --rm \
#     -v "$PWD/models:/models:ro" -v "$PWD/data:/data:ro" -v "$PWD/out:/out" \
#     ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml/mtb-resistotyper-ml-cli:vX.Y.Z \
#     predict --vcf /data/isolate.vcf --models /models --outdir /out
#
# `--models` may also point at the nested MODELS/<release>/<DRUG>/ layout
# deploy/webapp/bring-up.sh produces for more than one release; pass the
# specific release's subdirectory.

# --- build: resolve the pixi `prod` environment from pixi.toml/pixi.lock ---
FROM debian:bookworm-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*
# Pinned to what generated pixi.lock (lock format v7): an older pixi refuses
# to read it at all ("lock file version is 7, but only up to ... 6 is
# supported"), confirmed by trying 0.46.0 here first.
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.69.0 bash
ENV PATH="/root/.pixi/bin:${PATH}"

# Built and copied at the SAME absolute path in both stages, deliberately.
# pixi's installed console-script (`mtb-resistotyper-ml`) carries a shebang
# hardcoding this build-time path (`#!/opt/mtb/.pixi/envs/prod/bin/python3.12`)
# -- confirmed by building at /src and copying to /opt/pixi/prod first, which
# produced a script pointing at a python that no longer existed at runtime
# ("no such file or directory" on every invocation, including --version).
WORKDIR /opt/mtb
# Only what pixi needs to resolve/build the prod environment -- the wheel
# builds from this same copy of src/, and copying the whole repo would bust
# Docker's layer cache on every doc or test change.
COPY pixi.toml pixi.lock pyproject.toml README.md LICENSE ./
COPY src/ ./src/
# --locked: fail rather than silently re-resolve if pixi.lock and pixi.toml
# have drifted apart -- the same guarantee `--frozen` gives at the CLI, made
# the default so a stale lock cannot ship unnoticed.
RUN pixi install --environment prod --locked

# --- runtime: the resolved environment plus the reference and catalogue ---
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --create-home --uid 1000 mtb

# The resolved conda + PyPI environment pixi built -- bcftools, the package,
# gumpy/piezo/biopython==1.83 -- copied whole rather than re-resolved, so the
# runtime image never needs pixi or a package index itself. Same path as the
# build stage (see the WORKDIR comment above): the console-script shebang
# depends on it.
COPY --from=build /opt/mtb/.pixi/envs/prod /opt/mtb/.pixi/envs/prod
ENV PATH="/opt/mtb/.pixi/envs/prod/bin:${PATH}"

# Same source, same commit, as app/Dockerfile: the exact GenBank record the
# piezo/gumpy toolchain is built around, pinned to a commit rather than
# fetched from NCBI at build time.
ARG CATALOGUE_REF=master
ARG CATALOGUE_BASE=https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues
RUN mkdir -p /opt/reference /opt/catalogue \
 && curl -sSfL --retry 5 "${CATALOGUE_BASE}/${CATALOGUE_REF}/catalogues/NC_000962.3/NC_000962.3.gbk" \
      -o /opt/reference/NC_000962.3.gbk \
 && curl -sSfL --retry 5 "${CATALOGUE_BASE}/${CATALOGUE_REF}/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv" \
      -o /opt/catalogue/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
 && grep -q "^LOCUS" /opt/reference/NC_000962.3.gbk \
 && head -1 /opt/catalogue/*.csv | grep -q "GENBANK_REFERENCE"

ENV MTB_REFERENCE_GENBANK=/opt/reference/NC_000962.3.gbk \
    MTB_CATALOGUE=/opt/catalogue/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
    PYTHONUNBUFFERED=1
# No MTB_MODELS default: a bare `docker run` with no --models mount should
# say so, not silently score against nothing.

USER mtb
WORKDIR /home/mtb
ENTRYPOINT ["mtb-resistotyper-ml"]
CMD ["--help"]
