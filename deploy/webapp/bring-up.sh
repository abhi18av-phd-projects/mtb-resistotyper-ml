#!/usr/bin/env bash
# Build the serving image on aither and deploy it with the abc CLI.
#
# The image is built ON aither because "aither.local/..." is a naming
# convention for images already present in that host's docker daemon, not a
# registry that can be pushed to. Build elsewhere and the deploy places a job
# that can never pull.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MODELS_REPO="${MODELS_REPO:-$(cd "$REPO_ROOT/../mtb-resistotyper-ml-models" && pwd)}"
TAG="${TAG:-v0.1.0}"
HOST="${HOST:-sun-aither}"
CTX="${ABC_CLI_CONTEXT:-seedling-abhi}"

echo "==> building the runner wheel"
cd "$REPO_ROOT"
rm -rf dist && python3 -m build --wheel >/dev/null

echo "==> staging build context (runner + app + model bundles)"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/dist" "$STAGE/app"
cp dist/*.whl "$STAGE/dist/"
cp -R app/main.py app/vcf_to_garc.py app/static app/Dockerfile "$STAGE/app/"
# The bundles ship inside the image at a pinned models-repo revision. The
# artefact contract lets them move independently of the runner; pinning the
# revision here is what makes a served prediction attributable to one of them.
cp -R "$MODELS_REPO/models" "$STAGE/models"
git -C "$MODELS_REPO" rev-parse HEAD > "$STAGE/models/MODELS_REVISION"
mv "$STAGE/app/Dockerfile" "$STAGE/Dockerfile"

echo "==> shipping context to $HOST"
tar -C "$STAGE" -czf "$STAGE.tgz" .
ssh "$HOST" 'rm -rf ~/mtb-webapp-build && mkdir -p ~/mtb-webapp-build'
scp -q "$STAGE.tgz" "$HOST:~/mtb-webapp-build/ctx.tgz"
rm -f "$STAGE.tgz"

echo "==> building aither.local/mtb-resistotyper-webapp:$TAG on $HOST"
ssh "$HOST" "cd ~/mtb-webapp-build && tar xzf ctx.tgz && \
             docker build -t aither.local/mtb-resistotyper-webapp:$TAG ."

echo "==> deploying"
cd "$REPO_ROOT/deploy/webapp"
ABC_CLI_CONTEXT="$CTX" abc app validate
ABC_CLI_CONTEXT="$CTX" abc app deploy
ABC_CLI_CONTEXT="$CTX" abc app show mtb-resistotyper
