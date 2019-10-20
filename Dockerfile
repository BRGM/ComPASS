FROM registry.gitlab.inria.fr/charms/compass/run-environment:latest

ARG WHEEL_TAG

COPY wheel/*.whl /wheel/
RUN pip3 install /wheel/ComPASS-*${WHEEL_TAG}*.whl \
 && rm -rf /wheel

WORKDIR /localfs
