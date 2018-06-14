From debian:buster as builder
# as builder remove because gitlab runner version is older than this feature.
#to uncomment when building @BRGM
#ENV http_proxy='http://proxy.brgm.fr:8080/'
RUN  apt-get update && apt-get install --yes build-essential gcc gfortran cmake libmetis-dev python3 mpi-default-dev petsc-dev python3-dev
WORKDIR /build/
RUN mkdir /source/
COPY ./ /source/ComPASS-develop
#RUN cd /source && tar xzf ComPASS-develop.tar.gz && ls -al /source/

RUN CC=mpicc cmake ../source/ComPASS-develop && make



From debian:buster-slim
ENV PYTHONPATH=/build/meshtoolsModule:/build/compassModule
#to uncomment when building @BRGM
#ENV http_proxy='http://proxy.brgm.fr:8080/'
RUN apt-get update && apt-get install --yes --no-install-recommends python3-mpi4py python3-numpy metis libpetsc-real3.8 && apt-get clean 
#&& rm -rf /var/lib/apt/lists/*
WORKDIR /build/
COPY ./docker/script/docker_entrypoint.sh /
COPY --from=builder /source/ComPASS-develop/python ./compassModule
COPY --from=builder /source/ComPASS-develop/thirdparty/meshtools/python ./meshtoolsModule
VOLUME [/data]


#to uncomment when bash script will be done.
ENTRYPOINT ["/bin/bash","/docker_entrypoint.sh"]