FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest

RUN apt-get update && apt-get install --yes --no-install-recommends \
    cmake-curses-gui \
    python3-scipy \
    python3-matplotlib \
    ipython3 \
 && apt-get clean

USER compass
WORKDIR /localfs
