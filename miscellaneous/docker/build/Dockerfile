From debian:buster as builder
RUN  apt-get update && apt-get install --yes build-essential gcc gfortran cmake libmetis-dev libcgal-dev python3 mpi-default-dev petsc-dev python3-dev
WORKDIR /build/
ENTRYPOINT ["/bin/bash"]

