"""Read VCFs a researcher already uploaded, rather than asking for them twice.

The cluster has an upload path: `abc data upload` from the command line, and
tusd from a browser, both landing in the group bucket. A service that insists on
its own multipart form makes the researcher upload a second copy of a file the
platform already holds, and that copy is unattributed and unretained.

So the object store is the primary input. The direct upload stays as a fallback
for someone with a file and no cluster account, but the flow the platform is
built around is: upload once, then point at it.

Credentials arrive as AWS_* in the task environment, minted per-app by the
`data:` block in abc-app.yaml. Nothing is stored here and nothing is written
back: the service reads what it was pointed at, and refuses anything outside the
prefixes it was configured with.
"""

from __future__ import annotations

import os
from pathlib import Path

BUCKET = os.environ.get("MTB_BUCKET", "su-mbhg-bioinformatics")
# Where uploads land. `abc data upload` writes under user/<slot>/; the tusd
# mover writes user/<slot>/data/. Both are listed so a file is findable however
# it arrived.
PREFIXES = [p for p in os.environ.get(
    "MTB_PREFIXES", "user/abhinav/,users/abhi/,shared/").split(",") if p]
SUFFIXES = (".vcf", ".vcf.gz", ".bcf")


class StorageUnavailable(RuntimeError):
    pass


def _client():
    try:
        import boto3
        from botocore.config import Config
    except ImportError as exc:
        raise StorageUnavailable("boto3 is not installed") from exc
    endpoint = os.environ.get("ABC_MINIO_ENDPOINT") or os.environ.get("AWS_ENDPOINT_URL")
    if not endpoint:
        raise StorageUnavailable("no ABC_MINIO_ENDPOINT in the environment")
    if not os.environ.get("AWS_ACCESS_KEY_ID"):
        raise StorageUnavailable("no AWS credentials in the environment")
    return boto3.client("s3", endpoint_url=endpoint,
                        config=Config(s3={"addressing_style": "path"},
                                      retries={"max_attempts": 3}))


def _allowed(key: str) -> bool:
    """A key must sit under a configured prefix.

    Without this the service would fetch any object the app's credentials can
    reach on behalf of anyone who can reach the page, which turns a prediction
    endpoint into a general-purpose read oracle for the group bucket.
    """
    return any(key.startswith(p) for p in PREFIXES)


def list_vcfs(limit: int = 200) -> list[dict]:
    s3 = _client()
    out: list[dict] = []
    for prefix in PREFIXES:
        token = None
        while len(out) < limit:
            kw = {"Bucket": BUCKET, "Prefix": prefix, "MaxKeys": 1000}
            if token:
                kw["ContinuationToken"] = token
            page = s3.list_objects_v2(**kw)
            for o in page.get("Contents", []):
                if o["Key"].lower().endswith(SUFFIXES):
                    out.append({"key": o["Key"], "size": o["Size"],
                                "modified": o["LastModified"].isoformat()})
            token = page.get("NextContinuationToken")
            if not token:
                break
    out.sort(key=lambda r: r["modified"], reverse=True)
    return out[:limit]


def fetch(key: str, dest_dir: Path) -> Path:
    if not _allowed(key):
        raise StorageUnavailable(
            f"{key!r} is outside the prefixes this service may read "
            f"({', '.join(PREFIXES)})")
    s3 = _client()
    dest = dest_dir / Path(key).name
    dest.parent.mkdir(parents=True, exist_ok=True)
    s3.download_file(BUCKET, key, str(dest))
    return dest
