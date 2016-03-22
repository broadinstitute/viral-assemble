#!/bin/bash
set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
    brew update > /dev/null
    brew tap homebrew/versions
    brew unlink gcc
    brew unlink boost
    brew install gcc49 
    brew install boost 
    brew reinstall boost --with-python
    hash -r
fi