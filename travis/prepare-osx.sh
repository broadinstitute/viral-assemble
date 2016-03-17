#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    brew update
    #brew install python3

    # Install some custom requirements on OS X
       # e.g. brew install pyenv-virtualenv

        case "${TRAVIS_OSX_PYTHON_VERSION}" in
            py27)
                # Install some custom Python 2.7 requirements on OS X
                brew install python
                ;;
            py35)
                # Install some custom Python 3.5 requirements on OS X
                brew install python3
                ;;
        esac
fi