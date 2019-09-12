#!/bin/bash
#
# This script is intended to facilitate a development environment for working
# on viral-assemble. It is intended to be run within a viral-core container
# that has a git checkout of viral-assemble mounted into the container as
# /opt/viral-ngs/viral-assemble.
#
# It should be run once after spinning up a plain viral-core container.
# It will install conda dependencies and create symlinks for code modules.
# You do not need to run or source this afterwards as long as you docker commit
# or otherwise save the state of your container.
#
# Not intended for use in docker build contexts (see Dockerfile for that)
#
# Must have $INSTALL_PATH and $VIRAL_NGS_PATH defined (from viral-core docker)

VIRAL_ASSEMBLE_PATH=/opt/viral-ngs/viral-assemble

$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_ASSEMBLE_PATH/requirements-conda.txt

ln -s $VIRAL_ASSEMBLE_PATH/assemble $VIRAL_ASSEMBLE_PATH/assembly.py $VIRAL_NGS_PATH
