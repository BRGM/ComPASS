FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest AS base

RUN apt-get update && apt-get install --yes --no-install-recommends \
    wget \
    build-essential \
    git \
    clang-format \
    gcc \
    gfortran \
    cmake \
    libopenmpi-dev \
    libmetis-dev \
    libpetsc-real3.12-dev \
    python3-mpi4py \
 && apt-get clean \
 && python3 -m pip install \
    pre-commit setuptools_git_versioning

COPY su-exec/su-exec.c /tmp/su-exec.c
RUN gcc -O3 -Wall -Werror /tmp/su-exec.c -o /usr/local/sbin/su-exec \
 && rm /tmp/su-exec.c \
 && mkdir -p /build

# Build petsc4py with local version of petsc
# We need python wheel and petsc dev packages here
# Canonical PETSc location in ubuntu
ENV PETSC_DIR=/usr/lib/petsc
RUN mkdir -p /build \
 && mkdir -p /wheel \
 && cd /build \
 && wget https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-3.12.0.tar.gz \
 && tar xf petsc4py-3.12.0.tar.gz \
 && cd petsc4py-3.12.0 \
 && python3 setup.py bdist_wheel \
 && cp -v dist/*.whl /wheel/ \
 && cd / \
 && rm -rf build
VOLUME /wheel

# install pre-commit environment
WORKDIR /dummy
COPY .pre-commit-config.yaml .
RUN git init \
 && pre-commit install \
 && pre-commit run --all \
 && rm -rf /dummy

# This is used by cmake to find PETSc from the openmpi link flags
ENV CC=mpicc

RUN mkdir -p /localfs
VOLUME /localfs
WORKDIR /localfs
