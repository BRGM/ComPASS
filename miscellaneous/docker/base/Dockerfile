FROM registry.gitlab.inria.fr/charms/meshtools:master AS meshtools

FROM ubuntu:19.04

# We need python3-dev to have dynamic libraries at runtime
RUN apt-get update && apt-get install --yes --no-install-recommends \
    python3-dev \
    python3-numpy \
    python3-setuptools \
    python3-pytest-xdist \
    python3-pip \
 && apt-get clean

# The following is due to an output error with openmpi3 \
# && export OMPI_MCA_btl_vader_single_copy_mechanism=none

# Meshtools installation
WORKDIR /wheels
COPY --from=meshtools /wheels/MeshTools-0.0.1-cp37-cp37m-linux_x86_64.whl .
RUN pip3 install MeshTools-0.0.1-cp37-cp37m-linux_x86_64.whl \
 && rm -rf /wheels

RUN mkdir -p /localfs
VOLUME /localfs
WORKDIR /localfs

