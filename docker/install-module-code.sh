#!/bin/bash
#

set -e -o pipefail

echo "PATH:              ${PATH}"
echo "INSTALL_PATH:      ${INSTALL_PATH}"
echo "VIRAL_NGS_PATH:    ${VIRAL_NGS_PATH}"

cd $INSTALL_PATH/viral-assemble
ln -s assembly.py assemble $VIRAL_NGS_PATH