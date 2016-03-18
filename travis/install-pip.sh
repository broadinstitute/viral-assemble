#!/bin/bash
set -e

echo "CC: $CC" 
echo "CXX: $CXX" 
echo "which gxx: $(which gcc)" 
echo "gcc --version: $(gcc --version)"
env

echo "pip installing required python packages"
pip install -r requirements.txt

#PYVER=`python -V 2>&1 | cut -c 8`
PYVER=`echo $TRAVIS_PYTHON_VERSION | cut -c 1`
if [ "$PYVER" = "3" ] || [[ "$TRAVIS_OSX_PYTHON_VERSION" == py3* ]]; then
    echo "pip installing snakemake packages (py3 only)"
    pip install -r requirements-pipes.txt
fi

python --version

echo "pip installing test-related packages (coveralls, etc.)"
pip install -r requirements-tests.txt
