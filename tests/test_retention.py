"""Anonymous uploads must expire, and the expiry must be verifiable.

These run without a bucket: what is asserted is the policy, not the network.
A sweep that deletes the wrong thing, or a key scheme that lets one anonymous
uploader read another's file, are both silent failures in production.
"""

from __future__ import annotations

import datetime
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "app"))
import storage


def test_listing_is_off_unless_a_prefix_is_named(monkeypatch) -> None:
    """The dangerous default is an open listing; assert we do not have it."""
    monkeypatch.setattr(storage, "LIST_PREFIXES", [])
    with pytest.raises(storage.ListingDisabled):
        storage.list_vcfs()


def test_reads_are_confined_to_configured_prefixes(monkeypatch) -> None:
    monkeypatch.setattr(storage, "READ_PREFIXES", ["users/abhi/", "anon-uploads/"])
    assert storage._allowed("users/abhi/x.vcf")
    assert storage._allowed("anon-uploads/2026-09-04/tok/x.vcf")
    # The classic escapes: a sibling prefix, and a prefix that merely starts the
    # same way as an allowed one.
    assert not storage._allowed("user/someone-else/x.vcf")
    assert not storage._allowed("../users/abhi/x.vcf")
    assert not storage._allowed("secrets/keys.json")


def test_anonymous_keys_are_dated_and_unguessable(monkeypatch) -> None:
    """One anonymous uploader must not reach another's file by guessing a name."""
    captured = {}

    class FakeS3:
        def put_object(self, **kw):
            captured.update(kw)

    monkeypatch.setattr(storage, "_client", lambda: FakeS3())
    monkeypatch.setattr(storage, "ANON_PREFIX", "anon-uploads/")
    k1 = storage.store_anonymous(b"##fileformat=VCFv4.2\n", "sample.vcf")
    k2 = storage.store_anonymous(b"##fileformat=VCFv4.2\n", "sample.vcf")
    today = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")

    assert k1.startswith(f"anon-uploads/{today}/") and k1.endswith("/sample.vcf")
    assert k1 != k2, "identical filenames must not collide into one key"
    token = k1.split("/")[2]
    assert len(token) >= 12, "the token is what makes the key unguessable"
    assert captured["Metadata"]["uploaded"] == "anonymous"


def test_a_traversing_filename_cannot_escape_the_prefix(monkeypatch) -> None:
    stub = type("F", (), {"put_object": lambda *a, **k: None})
    monkeypatch.setattr(storage, "_client", lambda: stub())
    key = storage.store_anonymous(b"x", "../../etc/passwd")
    assert key.startswith("anon-uploads/") and ".." not in key


def test_oversized_uploads_are_refused(monkeypatch) -> None:
    monkeypatch.setattr(storage, "ANON_MAX_BYTES", 10)
    with pytest.raises(storage.StorageUnavailable):
        storage.store_anonymous(b"x" * 11, "big.vcf")


def test_sweep_deletes_only_what_is_past_retention(monkeypatch) -> None:
    now = datetime.datetime.now(datetime.timezone.utc)
    objects = [
        {"Key": "anon-uploads/old/a.vcf", "Size": 100,
         "LastModified": now - datetime.timedelta(hours=48)},
        {"Key": "anon-uploads/fresh/b.vcf", "Size": 200,
         "LastModified": now - datetime.timedelta(hours=1)},
    ]
    deleted = []

    class FakeS3:
        def list_objects_v2(self, **kw):
            assert kw["Prefix"] == storage.ANON_PREFIX, "the sweep must be prefix-scoped"
            return {"Contents": objects}

        def delete_objects(self, Bucket, Delete):
            deleted.extend(o["Key"] for o in Delete["Objects"])

    monkeypatch.setattr(storage, "_client", lambda: FakeS3())
    r = storage.sweep_anonymous(max_age_hours=24)
    assert deleted == ["anon-uploads/old/a.vcf"], "the fresh upload must survive"
    assert r["expired"] == 1 and r["deleted"] == 1 and r["scanned"] == 2

    deleted.clear()
    r = storage.sweep_anonymous(max_age_hours=24, dry_run=True)
    assert deleted == [] and r["expired"] == 1 and r["deleted"] == 0


def test_lifecycle_rule_preserves_other_rules(monkeypatch) -> None:
    """Applying our expiry must not silently drop somebody else's rule."""
    put = {}

    class FakeS3:
        def get_bucket_lifecycle_configuration(self, Bucket):
            return {"Rules": [{"ID": "someone-elses", "Status": "Enabled"}]}

        def put_bucket_lifecycle_configuration(self, Bucket, LifecycleConfiguration):
            put.update(LifecycleConfiguration)

    monkeypatch.setattr(storage, "_client", lambda: FakeS3())
    storage.ensure_lifecycle(days=1)
    ids = [r["ID"] for r in put["Rules"]]
    assert "someone-elses" in ids and storage.LIFECYCLE_RULE_ID in ids
    ours = next(r for r in put["Rules"] if r["ID"] == storage.LIFECYCLE_RULE_ID)
    assert ours["Filter"]["Prefix"] == storage.ANON_PREFIX
    assert ours["Expiration"]["Days"] == 1
