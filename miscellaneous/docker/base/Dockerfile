FROM ubuntu:20.04

# We need python3-dev to have dynamic libraries at runtime
RUN apt-get update && apt-get install --yes --no-install-recommends \
    python3-dev \
    python3-click \
    python3-numpy \
    python3-setuptools \
    python3-pytest-xdist \
    python3-pip \
    python3-wheel \
    python3-matplotlib \
 && apt-get clean

# The following is due to an output error with openmpi3 \
# && export OMPI_MCA_btl_vader_single_copy_mechanism=none

# Meshtools installation
ARG MESHTOOLS_WHEEL_PATH
WORKDIR /wheel
COPY $MESHTOOLS_WHEEL_PATH .
RUN pip3 install ${MESHTOOLS_WHEEL_PATH##*/} \
 && rm -rf /wheel

RUN mkdir -p /localfs
VOLUME /localfs
WORKDIR /localfs
