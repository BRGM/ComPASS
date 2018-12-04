FROM registry.gitlab.inria.fr/charms/meshtools:master AS meshtools

FROM debian:buster-slim

COPY ./set-compass-ids.sh /etc/local/sbin/set-compass-ids.sh
RUN chown root:root /etc/local/sbin/set-compass-ids.sh && \
    chmod 755 /etc/local/sbin/set-compass-ids.sh

COPY ./add-compass-user.sh /etc/local/sbin/add-compass-user.sh
RUN chown root:root /etc/local/sbin/add-compass-user.sh && \
    chmod 700 /etc/local/sbin/add-compass-user.sh && \
    /etc/local/sbin/add-compass-user.sh
RUN mkdir -p /localfs
VOLUME /localfs

COPY ./customize-bashrc.sh /etc/customize-bashrc.sh
RUN echo "\n# ComPASS specific" >> /etc/bash.bashrc && \
    cat /etc/local/sbin/set-compass-ids.sh >> /etc/bash.bashrc && \
    cat /etc/customize-bashrc.sh >> /etc/bash.bashrc && \
    rm /etc/customize-bashrc.sh

# We install python3.7-dev because python3-scipy python3-matplotlib 
# that are installed in work environment depend on python-3.7
# Is this an inconsistency in debian buster?
RUN apt-get update && apt-get install --yes --no-install-recommends \
    python3.7-dev \
    python3-numpy \
    python3-setuptools \
    wget \
 && apt-get clean 

# We need to use python3.7
# python3.6 is here because of unconsistency between debian buster packages
# integers are priorities we leave it 1 empty for later install of python2.7 coming with dependencies
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.6 2 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.7 3

# We install pip from sources as python3-pip package uses python3.6
RUN apt-get install --yes --no-install-recommends \
    wget \
 && apt-get clean \
 && wget --no-check-certificate https://files.pythonhosted.org/packages/45/ae/8a0ad77defb7cc903f09e551d88b443304a9bd6e6f124e75c0fbbf6de8f7/pip-18.1.tar.gz \
 && tar xf pip-18.1.tar.gz \
 && cd pip-18.1 \
 && python setup.py install \
 && cd .. \
 && rm -rf pip-18.1*

# Meshtools installation
WORKDIR /wheels
COPY --from=meshtools /wheels/MeshTools-0.0.1-cp37-cp37m-linux_x86_64.whl .
RUN pip install MeshTools-0.0.1-cp37-cp37m-linux_x86_64.whl

WORKDIR /localfs
RUN rm -rf /wheels
