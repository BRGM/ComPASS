FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest AS base

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    gcc \
    gfortran \
    cmake \
    libopenmpi-dev \
    libmetis-dev \
    libpetsc-real3.9-dev \
 && apt-get clean

COPY su-exec/su-exec.c /tmp/su-exec.c
RUN gcc -O3 -Wall -Werror /tmp/su-exec.c -o /usr/local/sbin/su-exec \
 && rm /tmp/su-exec.c \
 && mkdir -p /build
