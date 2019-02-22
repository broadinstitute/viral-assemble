#!/bin/bash

set -eu -o pipefail -x

# log in to DNAnexus
dx login --token "$DX_API_TOKEN" --noprojects
dx select $DX_PROJECT
dx pwd
dx ls -l
dx download file-FVQ74G80f5zf4V5vK0GjB1q5

pytest --cov-append $PYTEST_EXTRA_OPTS test/unit/test_workflow_utils.py::test_starting_cromwell_server

rc=$?; if [[ $rc != 0 ]]; then sleep 10; exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
