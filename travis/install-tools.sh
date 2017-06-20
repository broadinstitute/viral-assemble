#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # some conda packages dont exist on OSX
    cat requirements-conda.txt | grep -v kraken > $HOME/requirements-conda.txt
else
    # for linux, just use requirements-conda as-is
    cp requirements-conda.txt $HOME
fi

# Set to conda's java
export JAVA_HOME="$(pwd)/tools/conda-tools/default/jre"

echo "Installing and validating bioinformatic tools"
export CONDA_ENVS_PATH=tools/conda-cache:tools/conda-tools/default

for i in $(seq 3); do
  conda create -y -q -m -c broad-viral -c r -c bioconda -c conda-forge -c defaults -p tools/conda-tools/default --file $HOME/requirements-conda.txt python="$TRAVIS_PYTHON_VERSION" && break
  sleep 5
done

echo 'Sourcing default environment'
source activate tools/conda-tools/default
conda info -a # for debugging

./install_tools.py
