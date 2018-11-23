From debian:buster as builder
RUN  apt-get update && apt-get install --yes --no-install-recommends build-essential gcc gfortran cmake libmetis-dev libcgal-dev python3 mpi-default-dev petsc-dev python3-dev && apt-get clean
WORKDIR /build/
ENTRYPOINT ["/bin/bash"]

