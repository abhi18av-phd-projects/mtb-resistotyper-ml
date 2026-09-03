# The serving image. Deliberately not the training image.
#
# Scoring a deployed model is a scaler and a coefficient vector, so the runner
# itself needs nothing but the standard library. What makes this image large is
# VCF annotation: turning genomic coordinates into GARC gene@mutation names needs
# gumpy and an H37Rv GenBank. That cost buys names that mean something outside
# this repository, which a REF-POS-ALT coordinate string does not.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# The reference is baked in, not fetched at request time. A service that reaches
# out to NCBI to answer a prediction fails whenever NCBI does, and silently
# changes its answers whenever the record is revised.
RUN mkdir -p /opt/reference && curl -sSfL --retry 5 \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000962.3&rettype=gbwithparts&retmode=text" \
  -o /opt/reference/NC_000962.3.gbk \
 && test -s /opt/reference/NC_000962.3.gbk \
 && grep -q "^LOCUS" /opt/reference/NC_000962.3.gbk

RUN pip install --no-cache-dir \
      fastapi==0.115.6 "uvicorn[standard]==0.34.0" python-multipart==0.0.20 \
      gumpy==1.3.8 \
      # gumpy reads SeqFeature.strand, which Biopython deprecated in 1.81 and
      # removed in 1.85. An unpinned install therefore resolves to a Biopython
      # that parses the GenBank file and then fails on every gene, so the pin is
      # load-bearing rather than tidy.
      "biopython==1.83"

COPY dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -f /tmp/*.whl

WORKDIR /srv
COPY models /opt/models
COPY app/main.py app/vcf_to_garc.py /srv/
COPY app/static /srv/static
ENV MTB_MODELS=/opt/models \
    MTB_REFERENCE_GENBANK=/opt/reference/NC_000962.3.gbk \
    PYTHONUNBUFFERED=1

EXPOSE 8080
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]
