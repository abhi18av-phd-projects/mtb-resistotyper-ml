# Platform automations

Job definitions for the obligations that belong to the platform rather than to
the prediction service: data arriving, and data not staying.

They are `abc job run` scripts. The `#ABC` preamble carries the placement and
resources, and the body is what runs inside the container, so a definition is
readable as the command somebody would otherwise type.

| Job | Deletes? | What it is for |
| --- | --- | --- |
| `retention-verify.sh` | no | Is the retention promise being kept? Reports the bucket rule and what a sweep would remove. Exits non-zero when no rule is present. |
| `lifecycle-ensure.sh` | no | Applies the store-side expiry rule. Idempotent; run after every deployment. |
| `retention-sweep.sh` | **yes** | Deletes uploads past their window. The backstop for a rule that was never applied. |

```bash
abc job run abc-cluster/retention-verify.sh --wait --logs
```

## Why these are jobs and not application code

The service already runs a daily retention loop in-process. That was enough to
make the behaviour exist and not enough to make it a guarantee, for three
reasons a job definition fixes.

It was invisible. Retention ran as a side effect of the web process being up, so
"are uploads actually expiring" could only be answered by reading `/health` on a
service that happened to be running. As jobs it has a submission, an exit status
and a log, like any other work on the cluster.

It was coupled to unrelated uptime. A promise about somebody else's data should
not depend on whether a prediction endpoint is healthy. The store-side rule is
the real mechanism; `lifecycle-ensure` is what puts it there.

It could not be verified without being trusted. `retention-verify` reports
without deleting, so the guarantee can be checked by someone who is not willing
to run a destructive job, including in a review.

The in-app loop is left in place as a third fallback. It can be disabled once
these are scheduled, and it should be: two mechanisms that both silently do
nothing look exactly like two mechanisms that both work.

## The ordering that matters

`ensure` → `verify` → (only if verify shows a backlog) `sweep`.

The dangerous failure here is not a sweep that deletes too much. It is a rule
that was never applied and a sweep nobody ran, which leaves uploads promised a
24-hour life sitting in a bucket indefinitely, and which looks from the outside
exactly like a system that is working. `verify` exists to make that state
noisy: it is the only one of the three that fails.

## What this does not yet solve

**There is no periodic job type.** `abc job run` submits one-shot jobs; the
platform has no cron or schedule primitive, so recurrence needs an external
trigger. Until then these are run after deployment and on demand, and the
in-app loop remains the only thing running them daily. That is the honest
status: the automations are declared and auditable, not yet scheduled.

**The service manifest still pins the node-local tag.** These jobs now name the
registry image, so they place on any node. `deploy/webapp/abc-app.yaml` beside
them still names `aither.local/…`, which exists only in one node's image cache;
until that is changed too, the service and its automations are pinned to
different copies of the same build. They live in one directory so that
divergence is visible rather than discovered during an incident.

**Credentials come from the job's identity.** The scripts read the bucket from
the image's environment and expect the platform to supply object-store
credentials, the same way the service receives them. A job submitted by a user
without write access to the bucket will fail on `sweep`, which is the correct
outcome.

## Configuration

Inherited from the service image, overridable per submission:

| Variable | Meaning |
| --- | --- |
| `MTB_BUCKET` | bucket holding the uploads |
| `MTB_ANON_PREFIX` | prefix anonymous uploads land under (`anon-uploads/`) |
| `MTB_ANON_RETENTION_HOURS` | retention window, default 24 |
