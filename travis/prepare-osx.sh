#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    brew update > /dev/null
    brew tap homebrew/versions
    brew unlink gcc
    brew unlink boost
    brew unlink cmake
    brew install gcc49 
    brew install boost 
    brew install cmake 

    #export CC=gcc-4.8
    #export CXX=g++-4.8

    hash -r

    #pip install --upgrade virtualenv
    #pip install --upgrade -e git+https://github.com/pypa/virtualenv.git#egg=virtualenv

    # case "${TRAVIS_OSX_PYTHON_VERSION}" in
    #     py27)
    #         # Install some custom Python 2.7 requirements on OS X
    #         #brew install python
    #         if [ ! -e $HOME/virtualenv/bin/activate ]; then
    #             cd $HOME
    #             #virtualenv venv --no-setuptools --no-pip --no-wheel
    #         fi
    #         ;;
    #     py35)
    #         # Install some custom Python 3.5 requirements on OS X
    #         #brew install python3
    #         if [ ! -e $HOME/virtualenv/bin/activate ]; then
    #             cd $HOME
    #             #pyvenv virtualenv
    #             #virtualenv venv -p $(which python3) --no-setuptools --no-pip --no-wheel
    #         fi
    #         ;;
    # esac
fi