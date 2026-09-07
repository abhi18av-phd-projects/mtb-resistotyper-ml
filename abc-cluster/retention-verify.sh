#!/usr/bin/env bash
#ABC --name=mtb-retention-verify
#ABC --driver=docker
#ABC --driver.config.image=aither.local/mtb-resistotyper-webapp:v0.2.5
#ABC --cores=1
#ABC --mem=512M
#ABC --time=00:05:00
#
# Reports whether the retention promise is being kept. Deletes nothing, so this
# is the one to run often and the one to run first after any deployment.
#
# Exits non-zero when the bucket carries no lifecycle rule. That is the failure
# worth alarming on: the uploads are still there, but the guarantee that they
# expire is not, and nothing else makes that visible.
set -euo pipefail
cd /srv
python -m retention verify
