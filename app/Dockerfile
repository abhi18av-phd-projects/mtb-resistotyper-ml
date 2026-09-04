# The serving image. Deliberately not the training image.
#
# Scoring a deployed model is a scaler and a coefficient vector, so the runner
# needs nothing but the standard library. What makes this image large is
# everything around the scoring: bcftools to normalise a real caller's VCF,
# gumpy plus an H37Rv GenBank to turn coordinates into GARC names, and piezo
# plus the WHO catalogue for Layer 1.
FROM python:3.11-slim

# bcftools does the VCF normalisation. Hand-rolled parsing was tried and did not
# survive contact with real MAGMA output; every transformation this service
# applies is a documented bcftools operation.
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates bcftools tabix \
 && rm -rf /var/lib/apt/lists/*

# Reference and catalogue come from the Oxford catalogue repository, pinned to a
# commit, rather than from NCBI at build time. It is the exact GenBank record the
# piezo/gumpy toolchain is built around, and it does not change under us when
# NCBI revises an annotation or when the network is unavailable.
ARG CATALOGUE_REF=master
ARG CATALOGUE_BASE=https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues
RUN mkdir -p /opt/reference /opt/catalogue \
 && curl -sSfL --retry 5 "${CATALOGUE_BASE}/${CATALOGUE_REF}/catalogues/NC_000962.3/NC_000962.3.gbk" \
      -o /opt/reference/NC_000962.3.gbk \
 && curl -sSfL --retry 5 "${CATALOGUE_BASE}/${CATALOGUE_REF}/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv" \
      -o /opt/catalogue/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
 && grep -q "^LOCUS" /opt/reference/NC_000962.3.gbk \
 && head -1 /opt/catalogue/*.csv | grep -q "GENBANK_REFERENCE"

RUN pip install --no-cache-dir \
      fastapi==0.115.6 "uvicorn[standard]==0.34.0" python-multipart==0.0.20 \
      boto3==1.35.99 \
      gumpy==1.3.8 piezo==0.9.2 \
      # gumpy reads SeqFeature.strand, which Biopython deprecated in 1.81 and
      # removed in 1.85. An unpinned install resolves to a Biopython that parses
      # the GenBank file and then fails on every gene, so the pin is
      # load-bearing rather than tidy.
      "biopython==1.83"

COPY dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -f /tmp/*.whl

WORKDIR /srv
COPY models /opt/models
COPY app/main.py app/normalise.py app/storage.py app/vcf_to_garc.py /srv/
COPY app/static /srv/static
ENV MTB_MODELS=/opt/models \
    MTB_REFERENCE_GENBANK=/opt/reference/NC_000962.3.gbk \
    MTB_CATALOGUE=/opt/catalogue/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
    PYTHONUNBUFFERED=1

EXPOSE 8080
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]
