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
import re
from pathlib import Path

BUCKET = os.environ.get("MTB_BUCKET", "su-mbhg-bioinformatics")
# Two separate permissions, because they leak differently.
#
# READ_PREFIXES is what /predict may fetch. Knowing a key is the precondition,
# and keys are not guessable in practice.
#
# LIST_PREFIXES is what /objects may ENUMERATE, and enumeration is the dangerous
# one: object keys in a research group's bucket routinely carry sample and
# patient identifiers in the file name, so listing them discloses who was
# sequenced even when every object stays unread. On a publicly reachable
# deployment this must be a curated demo prefix or nothing at all. It defaults to
# empty: a listing endpoint that silently exposes the whole bucket because nobody
# set a variable is the wrong default.
READ_PREFIXES = [p for p in os.environ.get(
    "MTB_PREFIXES", "user/abhinav/,users/abhi/,shared/").split(",") if p]
LIST_PREFIXES = [p for p in os.environ.get("MTB_LIST_PREFIXES", "").split(",") if p]

# Kept as an alias so callers that only care about reading keep working.
PREFIXES = READ_PREFIXES
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
    return any(key.startswith(p) for p in READ_PREFIXES)


class ListingDisabled(RuntimeError):
    """No prefix is approved for enumeration on this deployment."""


def list_vcfs(limit: int = 200) -> list[dict]:
    if not LIST_PREFIXES:
        raise ListingDisabled(
            "object listing is not enabled on this deployment. Object keys can "
            "carry sample identifiers, so enumeration is opt-in per prefix "
            "(MTB_LIST_PREFIXES). Supply an object key directly instead.")
    s3 = _client()
    out: list[dict] = []
    for prefix in LIST_PREFIXES:
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


# ── Anonymous uploads ────────────────────────────────────────────────────────
#
# Anonymous upload is a deliberate testing affordance, and retention is the part
# that has to be got right rather than left to good intentions. An uploaded VCF
# is a clinical isolate: the file name alone routinely carries a sample or
# patient identifier, and nobody knows who uploaded it or under what basis, so
# there is no one to ask and no consent recorded. Keeping such a file
# indefinitely is the failure mode; keeping it briefly, under a rule that cannot
# quietly stop running, is the mitigation.
#
# Two mechanisms, deliberately. The bucket lifecycle rule is the primary one
# because it is declarative and server-side: it survives the app crashing, being
# redeployed, or being scaled to zero, none of which a cron job survives. The
# in-process sweeper is the backstop, because a lifecycle rule that was never
# applied, or was applied to the wrong prefix, fails silently and looks
# identical to one that is working.

ANON_PREFIX = os.environ.get("MTB_ANON_PREFIX", "anon-uploads/")
ANON_RETENTION_HOURS = int(os.environ.get("MTB_ANON_RETENTION_HOURS", "24"))
ANON_MAX_BYTES = int(os.environ.get("MTB_ANON_MAX_BYTES", str(64 * 1024 * 1024)))
LIFECYCLE_RULE_ID = "mtb-resistotyper-anon-expiry"


def anon_enabled() -> bool:
    return os.environ.get("MTB_ANON_UPLOADS", "").lower() in {"1", "true", "yes"}


def store_anonymous(data: bytes, filename: str) -> str:
    """Store one anonymous upload under a dated, unguessable key."""
    import datetime
    import secrets
    if len(data) > ANON_MAX_BYTES:
        raise StorageUnavailable(
            f"upload is {len(data)} bytes; the anonymous limit is {ANON_MAX_BYTES}")
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", Path(filename or "upload.vcf").name)[:120]
    day = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
    # The token is what stops one anonymous uploader reading another's file by
    # guessing a name. Listing of this prefix stays off for the same reason.
    key = f"{ANON_PREFIX}{day}/{secrets.token_urlsafe(12)}/{safe}"
    _client().put_object(Bucket=BUCKET, Key=key, Body=data,
                         Metadata={"uploaded": "anonymous",
                                   "expires-after-hours": str(ANON_RETENTION_HOURS)})
    return key


def ensure_lifecycle(days: int | None = None) -> dict:
    """Apply, and report, the bucket rule that expires anonymous uploads.

    Idempotent: it reads the existing configuration, replaces only this rule, and
    leaves any other rule on the bucket alone. Expiration is in whole days
    because S3 lifecycle has no finer unit; the sweeper enforces the tighter
    hour-level retention.
    """
    days = max(1, (days if days is not None else ANON_RETENTION_HOURS) // 24 or 1)
    s3 = _client()
    rules = []
    try:
        current = s3.get_bucket_lifecycle_configuration(Bucket=BUCKET)
        rules = [r for r in current.get("Rules", []) if r.get("ID") != LIFECYCLE_RULE_ID]
    except Exception:
        pass                      # no configuration yet, which is the common case
    rules.append({
        "ID": LIFECYCLE_RULE_ID,
        "Status": "Enabled",
        "Filter": {"Prefix": ANON_PREFIX},
        "Expiration": {"Days": days},
        "AbortIncompleteMultipartUpload": {"DaysAfterInitiation": 1},
    })
    s3.put_bucket_lifecycle_configuration(
        Bucket=BUCKET, LifecycleConfiguration={"Rules": rules})
    return {"rule": LIFECYCLE_RULE_ID, "prefix": ANON_PREFIX, "expire_days": days}


def lifecycle_state() -> dict:
    """Whether the expiry rule is actually on the bucket, not whether we asked."""
    try:
        cfg = _client().get_bucket_lifecycle_configuration(Bucket=BUCKET)
    except Exception as exc:
        return {"present": False, "error": str(exc)[:160]}
    for r in cfg.get("Rules", []):
        if r.get("ID") == LIFECYCLE_RULE_ID:
            return {"present": True, "status": r.get("Status"),
                    "prefix": r.get("Filter", {}).get("Prefix"),
                    "expire_days": r.get("Expiration", {}).get("Days")}
    return {"present": False}


def sweep_anonymous(max_age_hours: int | None = None, dry_run: bool = False) -> dict:
    """Delete anonymous uploads past their retention. The backstop, not the rule."""
    import datetime
    hours = max_age_hours if max_age_hours is not None else ANON_RETENTION_HOURS
    cutoff = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=hours)
    s3 = _client()
    doomed, scanned, freed = [], 0, 0
    token = None
    while True:
        kw = {"Bucket": BUCKET, "Prefix": ANON_PREFIX, "MaxKeys": 1000}
        if token:
            kw["ContinuationToken"] = token
        page = s3.list_objects_v2(**kw)
        for o in page.get("Contents", []):
            scanned += 1
            if o["LastModified"] < cutoff:
                doomed.append({"Key": o["Key"]})
                freed += o["Size"]
        token = page.get("NextContinuationToken")
        if not token:
            break
    deleted = 0
    if doomed and not dry_run:
        for i in range(0, len(doomed), 1000):
            s3.delete_objects(Bucket=BUCKET, Delete={"Objects": doomed[i:i + 1000]})
            deleted += len(doomed[i:i + 1000])
    return {"scanned": scanned, "expired": len(doomed), "deleted": deleted,
            "bytes_freed": freed, "cutoff": cutoff.isoformat(), "dry_run": dry_run}
