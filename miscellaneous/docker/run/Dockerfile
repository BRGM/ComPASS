FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest

ENV PYTHONPATH=/build/meshtoolsModule:/build/compassModule

RUN apt-get update && apt-get install --yes --no-install-recommends \
    python3-mpi4py \
    metis \
    libcgal-dev \
    libpetsc-real3.9 \
 && apt-get clean 

RUN mkdir -p /build
