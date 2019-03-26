#!/bin/bash
set -e

# the miniconda directory may exist if it has been restored from cache
if [ -d "$MINICONDA_DIR" ] && [ -e "$MINICONDA_DIR/bin/conda" ]; then
    echo "Miniconda install already present from cache: $MINICONDA_DIR"
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        # on OSX we need to rely on the conda Python rather than the Travis-supplied system Python
        # so conda has a higher precedence
        export PATH="$MINICONDA_DIR/bin:$PATH"
    else
        export PATH="$MINICONDA_DIR/bin:$PATH"
    fi
    hash -r
else # if it does not exist, we need to install miniconda
    rm -rf "$MINICONDA_DIR" # remove the directory in case we have an empty cached directory

    if [[ "$TRAVIS_PYTHON_VERSION" == 2* ]]; then
        if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
            curl -S https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh > miniconda.sh;
        else
            curl -S https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh > miniconda.sh;
        fi
     else
        if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
            curl -S https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh > miniconda.sh;
        else
            curl -S https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh;
        fi
    fi

    bash miniconda.sh -b -p "$MINICONDA_DIR"
    chown -R "$USER" "$MINICONDA_DIR"
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        # on OSX we need to rely on the conda Python rather than the Travis-supplied system Python
        # so conda has a higher precedence
        export PATH="$MINICONDA_DIR/bin:$PATH"
    else
        export PATH="$MINICONDA_DIR/bin:$PATH"
    fi
    hash -r
    conda update -y -c conda-canary conda # for pre-release conda
    EXTRA_PINS="openssl=1.1.1b"
    conda config --set always_yes yes --set changeps1 no --set remote_max_retries 6 #--set channel_priority strict
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --add channels broad-viral
    # Use recommendations from https://github.com/bioconda/bioconda-recipes/issues/13774
    conda update -y conda
    # conda config --set channel_priority strict
    conda install -y pycryptosat ${EXTRA_PINS}
    conda config --set sat_solver pycryptosat
    conda install --quiet -y openjdk ${EXTRA_PINS} # openjdk=8.0.152
fi

# update certs
conda update --quiet -y ca-certificates certifi #openssl pyopenssl
conda info -a # for debugging
