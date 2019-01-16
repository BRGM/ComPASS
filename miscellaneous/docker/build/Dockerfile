FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest AS base

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    gcc \
    gfortran \
    cmake \
    libmetis-dev \
    libcgal-dev \
    mpi-default-dev \
    libpetsc-real-dev \
 && apt-get clean

RUN mkdir -p /build
