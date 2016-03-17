#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    brew update
    brew install pyenv-virtualenv

    # Install some custom requirements on OS X
       # e.g. brew install pyenv-virtualenv

        case "${TRAVIS_OSX_PYTHON_VERSION}" in
            py27)
                # Install some custom Python 2.7 requirements on OS X
                brew install python
                if [ ! -e $HOME/virtualenv/bin/activate ]; then
                    cd $HOME
                    virtualenv virtualenv
                fi
                ;;
            py35)
                # Install some custom Python 3.5 requirements on OS X
                brew install python3
                if [ ! -e $HOME/virtualenv/bin/activate ]; then
                    cd $HOME
                    pyvenv-3.5 virtualenv
                fi
                ;;
        esac


fi

# if we're on OSX we have to activat the virtualenv ourselves
# at least until Travis supports Python builds nativelt
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    # if the virtual env exists but is not active
    if [ -e "$HOME/virtualenv/bin/activate" ] && [ -z "$VIRTUAL_ENV" ]; then
        # activate the virtualenv
        source $HOME/virtualenv/bin/activate
    fi
fi