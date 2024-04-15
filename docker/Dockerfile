ARG SDK_TAG="latest"
FROM registry.gitlab.inria.fr/compass/v4/compass-images/run:${SDK_TAG}

COPY dist/*.whl /wheel/
RUN python3 -m pip install --no-index -f /wheel ComPASS \
 && rm -rf /wheel
