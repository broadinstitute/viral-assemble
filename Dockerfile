FROM quay.io/broadinstitute/viral-core:2.1.10

LABEL maintainer "viral-ngs@broadinstitute.org"

ENV VIRAL_ASSEMBLE_PATH=$INSTALL_PATH/viral-assemble

COPY requirements-conda.txt $VIRAL_ASSEMBLE_PATH/
RUN $VIRAL_NGS_PATH/docker/install-conda-dependencies.sh $VIRAL_ASSEMBLE_PATH/requirements-conda.txt

# Copy all source code into the base repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_ASSEMBLE_PATH

# Link key bits of python code into the path
RUN ln -s $VIRAL_ASSEMBLE_PATH/assembly.py $VIRAL_ASSEMBLE_PATH/assemble $VIRAL_NGS_PATH

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; assembly.py --version"

CMD ["/bin/bash"]
