#!/bin/bash
#set -e

py.test -v -rsxX -n 2 --durations=50 --fixture-durations=20 --junit-xml=pytest.xml --testmon test/unit test/integration
cp .coverage.cumulative .coverage

