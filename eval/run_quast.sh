#!/usr/bin/env bash

set -e -o pipefail
source activate py27
/idi/sabeti-scratch/ilya/sw/quast-3.2/quast.py "$@"


