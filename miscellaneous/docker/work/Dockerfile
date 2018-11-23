From registry.gitlab.inria.fr/charms/meshtools:master as meshtools

From debian:buster as builder

# We install python3.7-dev because python3-scipy python3-matplotlib depends on python-3.7
# Is this an inconsistency in debian buster?
RUN  apt-get update && apt-get install --yes --no-install-recommends build-essential gcc gfortran cmake-curses-gui libmetis-dev libcgal-dev mpi-default-dev petsc-dev python3.7-dev python3-pip python3-mpi4py python3-numpy python3-scipy python3-matplotlib ipython3 && apt-get clean

WORKDIR /wheels/
COPY --from=meshtools /wheels/MeshTools-0.0.1-cp36-cp36m-linux_x86_64.whl .
RUN pip3 install MeshTools-0.0.1-cp36-cp36m-linux_x86_64.whl
RUN rm -rf /build /wheels

RUN mkdir /localfs
VOLUME [/localfs]
WORKDIR /localfs

RUN useradd --create-home -s /bin/bash compass && chown compass:compass /localfs
USER compass
RUN echo "export CC=mpicc" >> /home/compass/.bashrc

ENTRYPOINT ["/bin/bash"]

