#!/usr/bin/env bash
#ABC --name=mtb-retention-sweep
#ABC --driver=docker
#ABC --driver.config.image=aither.local/mtb-resistotyper-webapp:v0.2.5
#ABC --cores=1
#ABC --mem=512M
#ABC --time=00:15:00
#
# Deletes anonymous uploads past their retention window. DESTRUCTIVE.
#
# Run lifecycle-ensure first: if the store-side rule is in place this job
# normally finds nothing, and finding a backlog here means the rule was missing.
# Prefer retention-verify to see what this would remove before running it.
set -euo pipefail
cd /srv
python -m retention sweep
