FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest AS base

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    gcc \
    gfortran \
    cmake \
    libmetis-dev \
    libcgal-dev \
    mpi-default-dev \
    libpetsc-real3.9-dev \
 && apt-get clean

# We need to use python3.7
# python3.6 is here because of unconsistency between debian buster packages
# python2.7 comes along with dependencies
RUN update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.6 2 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.7 3

RUN mkdir -p /build
