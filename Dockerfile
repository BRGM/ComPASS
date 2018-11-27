FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest AS build

RUN mkdir -p /source
COPY ./ /source/ComPASS
WORKDIR /build
RUN CC=mpicc cmake ../source/ComPASS && make

FROM registry.gitlab.inria.fr/charms/compass/run-environment:latest

WORKDIR /build
RUN mkdir -p /build/ComPASS
COPY --from=build /source/ComPASS/python ./ComPASS/python
ENV PYTHONPATH=/build/ComPASS/python

WORKDIR /localfs
USER compass

#to uncomment when bash script will be done.
#ENTRYPOINT ["/bin/bash","/docker_entrypoint.sh"]
#ENTRYPOINT ["mpirun", "-n", "`nproc`", "python3.7"]
ENTRYPOINT ["python"]
