#From debian:buster as builder
#WORKDIR /build/
#RUN  apt-get update && apt-get install --yes build-essential gcc gfortran cmake libmetis-dev python3 mpi-default-dev petsc-dev python3-dev
#RUN apt-get install --yes libcgal-dev
#RUN CC=mpicc cmake ../source/ComPASS-develop && make -j `nproc`
#RUN make install
#RUN apt-get install --yes python3-mpi4py python3-numpy && apt-get clean

# check docker version if label builder is not accepted
From debian:buster as builder
RUN  apt-get update && apt-get install --yes build-essential gcc gfortran cmake libmetis-dev python3 mpi-default-dev petsc-dev python3-dev
WORKDIR /build/
RUN mkdir /source/
COPY ./ /source/ComPASS-develop
RUN CC=mpicc cmake /source/ComPASS-develop && make -j `nproc` install

From debian:buster-slim
ENV PYTHONPATH=/build/meshtoolsModule:/build/compassModule
RUN apt-get update && apt-get install --yes --no-install-recommends libpython3-dev python3-mpi4py python3-numpy metis libpetsc-real3.8 && apt-get clean 
#&& rm -rf /var/lib/apt/lists/*
WORKDIR /build/
#COPY ./docker/script/docker_entrypoint.sh /
COPY --from=builder /source/ComPASS-develop/python ./compassModule
COPY --from=builder /source/ComPASS-develop/thirdparty/meshtools/python ./meshtoolsModule

#VOLUME [/data]
#to uncomment when bash script will be done.
#ENTRYPOINT ["/bin/bash","/docker_entrypoint.sh"]

#RUN useradd --create-home -s /bin/bash compass
#USER compass
VOLUME ["/data"]

WORKDIR /data/
#ENTRYPOINT ["mpirun -n `nproc` python3"]
ENTRYPOINT ["python3"]
CMD [""]
