FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest

RUN apt-get update && apt-get install --yes --no-install-recommends \
    vim  \
    cmake-curses-gui \
 && apt-get clean

# Scorep compilation and installation
RUN mkdir -p /build/scorep \
 && cd /build/scorep \
 && wget --no-check-certificate https://www.vi-hps.org/cms/upload/packages/scorep/scorep-6.0.tar.gz \
 && tar xf scorep-6.0.tar.gz \
 && mkdir -p scorep-6.0/build \
 && cd scorep-6.0/build \
 && ../configure --prefix=/usr/local --with-mpi=openmpi --enable-shared \
 && make -j `nproc` \
 && make install \
 && make clean \
 && cd /build \
 && rm -rf scorep

# petsc4py installation
RUN pip3 install /wheel/petsc4py-*.whl

WORKDIR /localfs

COPY entrypoint.sh /etc/compass/entrypoint.sh
COPY entrypoint.py /etc/compass/entrypoint.py
ENTRYPOINT ["/etc/compass/entrypoint.sh"]
CMD [""]
