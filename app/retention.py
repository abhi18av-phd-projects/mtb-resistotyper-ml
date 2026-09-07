"""Retention automation, as a command rather than a side effect of serving.

The expiry of anonymous uploads is a platform obligation, not a feature of the
prediction service: it must hold whether or not that service is running, and it
must be auditable by someone who never opens the app. Running it inside the web
process made it invisible and coupled it to an unrelated uptime.

Three verbs, deliberately separate because they carry different risk:

    ensure   apply the bucket lifecycle rule. Idempotent, server-side, and the
             primary mechanism: the store keeps expiring objects even when
             nothing of ours is running.
    verify   report only. Says whether the rule is present and what a sweep
             WOULD delete. Deletes nothing, so it is safe to run anywhere and
             is the one to schedule most often.
    sweep    delete what is already past its retention. The backstop, for the
             case the rule was never applied -- which fails silently and looks
             exactly like a rule that is working.

`verify` exists because the dangerous failure here is not a sweep that deletes
too much. It is a rule that was never applied and a sweep nobody ran, leaving
uploads that were promised a 24-hour life sitting in a bucket indefinitely.
"""

from __future__ import annotations

import argparse
import json
import sys

import storage


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("verb", choices=["ensure", "verify", "sweep"])
    p.add_argument("--max-age-hours", type=int, default=None,
                   help="override the retention window for this run")
    p.add_argument("--expire-days", type=int, default=None,
                   help="lifecycle rule window, for `ensure`")
    a = p.parse_args()

    if not storage.anon_enabled():
        print(json.dumps({"verb": a.verb, "skipped": "anonymous uploads are disabled"}))
        return 0

    if a.verb == "ensure":
        out = storage.ensure_lifecycle(days=a.expire_days)
    elif a.verb == "verify":
        # A dry sweep plus the rule's state: together they answer "is the promise
        # being kept", which neither answers alone.
        out = {"lifecycle": storage.lifecycle_state(),
               "would_delete": storage.sweep_anonymous(max_age_hours=a.max_age_hours,
                                                       dry_run=True)}
        rule = (out["lifecycle"] or {}).get("expire_days")
        if rule is None:
            # Non-zero so a scheduler surfaces it. The objects are not yet lost;
            # the guarantee is.
            print(json.dumps(out, indent=2))
            print("FAIL: no lifecycle rule on the bucket — retention is unenforced",
                  file=sys.stderr)
            return 2
    else:
        out = storage.sweep_anonymous(max_age_hours=a.max_age_hours)

    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
