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
# mtb-resistotyper-ml-models (unsuffixed, MODELS_REPO below) is the development
# checkout: fast to iterate on, and the only set with no fe_version/fe_id
# recorded. It is for testing this script and the release selector, not for
# serving. MODELS_REPOS defaults instead to the two editions meant to reach
# Zenodo, so a bare run of this script ships what the manuscript and the
# deposit both point at; set MODELS_REPOS=$MODELS_REPO to serve the dev set.
MODELS_REPO="${MODELS_REPO:-$(cd "$REPO_ROOT/../mtb-resistotyper-ml-models" && pwd)}"
_RELEASE_SETS="$(cd "$REPO_ROOT/../mtb-resistotyper-ml-release-sets" && pwd)"
DEFAULT_MODELS_REPOS="$_RELEASE_SETS/mtb-resistotyper-ml-models-v3.4.0+fe1.0.3:$_RELEASE_SETS/mtb-resistotyper-ml-models-v2.1.2+fe1.0.3"
TAG="${TAG:-v0.2.10}"
REGISTRY="${REGISTRY:-ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml}"
IMAGE="mtb-resistotyper-webapp"
HOST="${HOST:-sun-aither}"
CTX="${ABC_CLI_CONTEXT:-seedling-abhi}"

echo "==> building the runner wheel"
cd "$REPO_ROOT"
# The interpreter is a parameter, not whatever `python3` happens to resolve to:
# on one shell it was a conda 3.11 with `build`, on another the system 3.9
# without it, and the script died on this line with its output sent to
# /dev/null -- so the deploy stopped before doing anything while a filtered
# log looked like it had merely been quiet.
PYTHON="${PYTHON:-python3}"
if ! "$PYTHON" -c 'import build' 2>/dev/null; then
  echo "error: $PYTHON has no 'build' module. Set PYTHON to an interpreter that does" >&2
  echo "       (python >= 3.10 with 'pip install build'), e.g. PYTHON=/path/to/python3.11" >&2
  exit 1
fi
rm -rf dist && "$PYTHON" -m build --wheel >/dev/null

echo "==> staging build context (runner + app + model bundles)"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/dist" "$STAGE/app"
cp dist/*.whl "$STAGE/dist/"
cp -R app/*.py app/static app/Dockerfile "$STAGE/app/"
# The bundles ship inside the image at a pinned models-repo revision. The
# artefact contract lets them move independently of the runner; pinning the
# revision here is what makes a served prediction attributable to one of them.
#
# MODELS_REPOS takes a colon-separated list. One entry keeps the flat
# MODELS/<DRUG>/ layout every deployment has had. Two or more nest as
# MODELS/<release>/<DRUG>/, which is the layout the release selector needs
# before it can offer a choice: with a single set the UI has one option and no
# comparison to draw. The directory is named from the bundles' own
# data.release, not from the repo's path, so a set stays identifiable after
# somebody renames a checkout.
MODELS_REPOS="${MODELS_REPOS:-$DEFAULT_MODELS_REPOS}"
IFS=':' read -r -a _repos <<< "$MODELS_REPOS"
if [ "${#_repos[@]}" -eq 1 ]; then
  cp -R "${_repos[0]}/models" "$STAGE/models"
  git -C "${_repos[0]}" rev-parse HEAD > "$STAGE/models/MODELS_REVISION"
else
  mkdir -p "$STAGE/models"
  for repo in "${_repos[@]}"; do
    rel="$(python3 -c "
import json,sys,glob
cards=sorted(glob.glob(sys.argv[1] + '/models/*/model_card.json'))
for c in cards:
    d=json.load(open(c))
    r=(d.get('data') or {}).get('release')
    if r:
        print(r); break
else:
    raise SystemExit('no data.release in any model_card under ' + sys.argv[1])
" "$repo")"
    echo "    $repo -> models/$rel"
    cp -R "$repo/models" "$STAGE/models/$rel"
    # A git checkout pins the set by commit. A directory of freshly trained
    # bundles is not a checkout, so pin it by content instead: a sha256 over
    # every file, which identifies the set as exactly as a revision does and
    # does not pretend to a commit that never existed.
    # A checkout whose models/ differs from HEAD is served as it is on disk, so
    # a bare commit would misstate what is deployed; record the commit AND the
    # content hash, marked dirty.
    if git -C "$repo" rev-parse HEAD >/dev/null 2>&1 && git -C "$repo" diff --quiet HEAD -- models; then
      git -C "$repo" rev-parse HEAD > "$STAGE/models/$rel/MODELS_REVISION"
    elif git -C "$repo" rev-parse HEAD >/dev/null 2>&1; then
      h=$( ( cd "$repo/models" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 ) | shasum -a 256 | awk '{print $1}')
      echo "$(git -C "$repo" rev-parse HEAD)+dirty sha256:$h" > "$STAGE/models/$rel/MODELS_REVISION"
    else
      ( cd "$repo/models" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256 ) \
        | shasum -a 256 | awk '{print "sha256:" $1}' > "$STAGE/models/$rel/MODELS_REVISION"
    fi
  done
fi
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
