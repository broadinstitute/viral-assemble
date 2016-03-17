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

echo "TRAVIS_BRANCH: $TRAVIS_BRANCH"
echo "TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"

if [ $TRAVIS_PULL_REQUEST != "false" -o $TRAVIS_BRANCH = "master" -o -n "$TRAVIS_TAG" ]; then
    echo "This is on master or is a pull request: executing long running tests..."
    nosetests -v \
        --logging-clear-handlers \
        --with-timer \
        --with-xunit --with-coverage \
        --cover-inclusive --cover-branches --cover-tests \
        --cover-package broad_utils,illumina,assembly,interhost,intrahost,ncbi,read_utils,reports,taxon_filter,tools,util \
        -w test/integration/

else
    echo "This is not a pull request: skipping long running tests."
fi
