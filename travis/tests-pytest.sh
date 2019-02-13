#!/bin/bash

set -eu -o pipefail

# log in to DNAnexus
dx login --token "$DX_API_TOKEN" --noprojects
dx select $DX_PROJECT

pytest --cov-append $PYTEST_EXTRA_OPTS test

rc=$?; if [[ $rc != 0 ]]; then sleep 10; exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
