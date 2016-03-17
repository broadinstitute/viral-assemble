#!/bin/bash
set -e

# if we're on OSX we have to activat the virtualenv ourselves
# at least until Travis supports Python builds nativelt
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    # if the virtual env exists but is not active
    if [ -e "$HOME/virtualenv/bin/activate" ] && [ -z "$VIRTUAL_ENV" ]; then
        # activate the virtualenv
        source $HOME/virtualenv/bin/activate
    fi
fi

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
