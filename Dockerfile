FROM quay.io/broadinstitute/viral-core:2.0.1

LABEL maintainer "viral-ngs@broadinstitute.org"

# to build:
#   docker build .
#
# to run:
#   docker run --rm <image_ID> "<command>.py subcommand"
#
# to run interactively:
#   docker run --rm -it <image_ID>
#
# to run with GATK and/or Novoalign:
#   Download licensed copies of GATK and Novoalign to the host machine (for Linux-64)
#   export GATK_PATH=/path/to/gatk/
#   export NOVOALIGN_PATH=/path/to/novoalign/
#   docker run --rm -v $GATK_PATH:/gatk -v $NOVOALIGN_PATH:/novoalign -v /path/to/dir/on/host:/user-data <image_ID> "<command>.py subcommand"

COPY requirements-conda.txt $VIRAL_NGS_PATH/
RUN $VIRAL_NGS_PATH/docker/install-viral-ngs.sh minimal

# Copy all of the source code into the repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $VIRAL_NGS_PATH/

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
RUN /bin/bash -c "set -e; echo -n 'viral-ngs version: '; assembly.py --version"

CMD ["/bin/bash"]
