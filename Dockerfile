# check docker version if label builder is not accepted
From registry.gitlab.inria.fr/charms/compass/build-environment:latest as builder
WORKDIR /build/
RUN mkdir /source/
COPY ./ /source/ComPASS-develop
#RUN cd /source && tar xzf ComPASS-develop.tar.gz && ls -al /source/

RUN CC=mpicc cmake ../source/ComPASS-develop && make

From registry.gitlab.inria.fr/charms/meshtools:master as meshtools

From registry.gitlab.inria.fr/charms/compass/run-environment:latest

WORKDIR /wheels/
COPY --from=meshtools /wheels/MeshTools-0.0.1-cp36-cp36m-linux_x86_64.whl .
RUN pip3 install MeshTools-0.0.1-cp36-cp36m-linux_x86_64.whl

ENV PYTHONPATH=/build/compassModule
WORKDIR /build/
#COPY ./docker/script/docker_entrypoint.sh /
COPY --from=builder /source/ComPASS-develop/python ./compassModule
VOLUME [/data]

WORKDIR /data

RUN useradd --create-home -s /bin/bash compass && chown compass:compass /data
USER compass

#to uncomment when bash script will be done.
#ENTRYPOINT ["/bin/bash","/docker_entrypoint.sh"]
ENTRYPOINT ["python3"]
#ENTRYPOINT ["mpirun", "-n", "`nproc`", "python3"]


