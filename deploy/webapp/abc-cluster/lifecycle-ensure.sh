#!/usr/bin/env bash
#ABC --name=mtb-lifecycle-ensure
#ABC --driver=docker
#ABC --driver.config.image=ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml/mtb-resistotyper-webapp:v0.2.5
#ABC --cores=1
#ABC --mem=512M
#ABC --time=00:05:00
#
# Applies the bucket rule that expires anonymous uploads. Idempotent, so it is
# safe to re-run and safe to run after every deployment.
#
# This is the PRIMARY mechanism: it is enforced by the object store, so uploads
# keep expiring even when the prediction service is stopped, redeployed or
# crash-looping. The sweep is only a backstop for the case this was never run.
set -euo pipefail
cd /srv
python -m retention ensure
