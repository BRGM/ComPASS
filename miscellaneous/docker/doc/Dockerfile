FROM registry.gitlab.inria.fr/charms/meshtools/run-environment:latest

# The following two lines are necessary for the configuration of tzdata
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update \
 && apt-get install --yes \
    python3-sphinx \
    python3-sphinx-rtd-theme \
    pandoc \
    texlive \
    doxygen \
    graphviz \
 && apt-get clean \
 && pip3 install --upgrade recommonmark

VOLUME /localfs
WORKDIR /localfs
CMD ["/bin/bash"]
