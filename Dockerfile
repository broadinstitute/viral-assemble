FROM quay.io/broadinstitute/viral-core:2.0.1

LABEL maintainer "viral-ngs@broadinstitute.org"

COPY requirements-conda.txt $INSTALL_PATH/viral-assemble/requirements-conda.txt
RUN $INSTALL_PATH/viral-assemble/docker/install-conda-dependencies.sh $INSTALL_PATH/viral-assemble/requirements-conda.txt

# Copy key bits of source code into the base repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $INSTALL_PATH/viral-assemble
RUN $INSTALL_PATH/viral-assemble/docker/install-module-code.sh

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; assembly.py --version"

CMD ["/bin/bash"]
