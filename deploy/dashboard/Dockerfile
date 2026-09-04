# The Layer-2 evaluation dashboard.
#
# The evaluation is NOT computed here. Scoring 7,767 isolates against eight model
# bundles is a batch job whose answer does not change between page loads, and a
# dashboard that recomputes it per request is slow, and worse, is a different
# number every time the models move underneath it. The results are baked in with
# the code that produced them, so a screenshot of this page can be traced to one
# cohort build and one scoring run.
FROM python:3.11-slim

RUN pip install --no-cache-dir "streamlit==1.41.1" "pandas==2.2.3"

WORKDIR /srv
COPY app.py /srv/app.py
COPY data /data
ENV MTB_EVAL_DATA=/data \
    PYTHONUNBUFFERED=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

EXPOSE 8501
# --server.baseUrlPath is set at run time from the app descriptor: the private and
# shared planes serve under /apps/<name>/, and Streamlit's SPA fetches its own
# assets by absolute path, so without the prefix every asset 404s behind the
# reverse proxy while the page itself loads.
CMD ["sh", "-c", "streamlit run app.py --server.port=8501 --server.address=0.0.0.0 ${STREAMLIT_BASE_PATH:+--server.baseUrlPath=$STREAMLIT_BASE_PATH}"]
