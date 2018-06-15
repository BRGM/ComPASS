
# check docker version if label builder is not accepted
From debian:buster as builder
RUN  apt-get update && apt-get install --yes build-essential gcc gfortran cmake libmetis-dev python3 mpi-default-dev petsc-dev python3-dev
WORKDIR /build/
RUN mkdir /source/
COPY ./ /source/ComPASS-develop
#RUN cd /source && tar xzf ComPASS-develop.tar.gz && ls -al /source/

RUN CC=mpicc cmake ../source/ComPASS-develop && make



From debian:buster-slim
ENV PYTHONPATH=/build/meshtoolsModule:/build/compassModule
RUN apt-get update && apt-get install --yes --no-install-recommends libpython3-dev python3-mpi4py python3-numpy metis libpetsc-real3.8 && apt-get clean 
WORKDIR /build/
COPY ./docker/script/docker_entrypoint.sh /
COPY --from=builder /source/ComPASS-develop/python ./compassModule
COPY --from=builder /source/ComPASS-develop/thirdparty/meshtools/python ./meshtoolsModule
VOLUME [/data]

WORKDIR /data

RUN useradd --create-home -s /bin/bash compass && chown compass:compass /data
USER compass

#to uncomment when bash script will be done.
#ENTRYPOINT ["/bin/bash","/docker_entrypoint.sh"]
#ENTRYPOINT ["python3"]
ENTRYPOINT ["mpirun", "-n", "`nproc`", "python3"]


