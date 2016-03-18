#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    brew update > /dev/null
    #brew install homebrew/versions/gcc49

    hash -r

    case "${TRAVIS_OSX_PYTHON_VERSION}" in
        py27)
            # Install some custom Python 2.7 requirements on OS X
            #brew install python
            if [ ! -e $HOME/virtualenv/bin/activate ]; then
                cd $HOME
                virtualenv venv --no-setuptools --no-pip --no-wheel
            fi
            ;;
        py35)
            # Install some custom Python 3.5 requirements on OS X
            #brew install python3
            if [ ! -e $HOME/virtualenv/bin/activate ]; then
                cd $HOME
                #pyvenv virtualenv
                virtualenv venv -p $(which python3) --no-setuptools --no-pip --no-wheel
            fi
            ;;
    esac
fi