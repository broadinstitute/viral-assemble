#!/bin/bash

#pytest --cov-append test/unit
pytest --log-cli-level=DEBUG -vvvsx test/unit/test_workflow_utils.py::test_git_annex_basic

rc=$?; if [[ $rc != 0 ]]; then sleep 10; exit $rc; fi
# sleep to allow logs to be printed without truncation in the event of error
