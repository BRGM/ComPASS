FROM registry.gitlab.inria.fr/charms/compass/run-environment:latest

COPY wheel/*.whl /wheel/
RUN pip3 install --upgrade --find-links /wheel/ --no-index ComPASS \
 && rm -rf /wheel

WORKDIR /localfs
