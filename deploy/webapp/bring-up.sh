#!/usr/bin/env bash
# Build the serving image on aither and deploy it with the abc CLI.
#
# The image is built on aither because the build context and the models are
# already there, and it is then PUSHED to GHCR. The push is what the manifest
# depends on: it previously named "aither.local/...", a convention for images
# already present in that one host's docker daemon, so the deploy worked only
# where the build happened to sit and would place an unpullable job anywhere
# else. Keep both tags -- the local one is what the build produces, the registry
# one is what anything else can actually pull.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MODELS_REPO="${MODELS_REPO:-$(cd "$REPO_ROOT/../mtb-resistotyper-ml-models" && pwd)}"
TAG="${TAG:-v0.2.5}"
REGISTRY="${REGISTRY:-ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml}"
IMAGE="mtb-resistotyper-webapp"
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
cp -R app/*.py app/static app/Dockerfile "$STAGE/app/"
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

echo "==> building $REGISTRY/$IMAGE:$TAG on $HOST"
ssh "$HOST" "cd ~/mtb-webapp-build && tar xzf ctx.tgz && \
             docker build -t aither.local/$IMAGE:$TAG -t $REGISTRY/$IMAGE:$TAG ."

# Pushed before the deploy, not after: abc app deploy places a job that pulls
# this tag, so a deploy that runs first can only succeed on the node holding
# the build -- the failure this push exists to remove.
echo "==> pushing $REGISTRY/$IMAGE:$TAG"
ssh "$HOST" "docker push -q $REGISTRY/$IMAGE:$TAG"

echo "==> deploying"
cd "$REPO_ROOT/deploy/webapp"
ABC_CLI_CONTEXT="$CTX" abc app validate
# --node-pool is not optional here. This context supplies no head_pool, so the
# generated job carries no node_pool and lands in "default", which holds no
# nodes: the deploy is accepted and then never places, reporting "No nodes were
# eligible for evaluation". The image also exists only in aither's local daemon,
# and aither is the platform pool, so the placement and the image agree.
ABC_CLI_CONTEXT="$CTX" abc app deploy --node-pool platform --health-timeout 4m
ABC_CLI_CONTEXT="$CTX" abc app show mtb-resistotyper
