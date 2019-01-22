FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest AS build

RUN mkdir -p /source
COPY ./ /source/ComPASS
WORKDIR /build
RUN CC=mpicc cmake ../source/ComPASS && make -j `nproc`

FROM registry.gitlab.inria.fr/charms/compass/run-environment:latest

RUN mkdir -p /build/ComPASS
COPY --from=build /source/ComPASS/python /build/ComPASS/python
ENV PYTHONPATH=/build/ComPASS/python

WORKDIR /localfs
