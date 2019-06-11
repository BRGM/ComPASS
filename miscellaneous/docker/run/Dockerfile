FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest AS build

FROM registry.gitlab.inria.fr/charms/compass/base-environment:latest

COPY --from=build /usr/local/sbin/su-exec /usr/local/sbin/su-exec

RUN apt-get update && apt-get install --yes --no-install-recommends \
    python3-mpi4py \
    metis \
    openmpi-bin \
    libpetsc-real3.10 \
    libcgal13 \
    python3-scipy \
 && apt-get clean 

RUN mkdir -p /build

COPY entrypoint.sh /etc/compass/entrypoint.sh
COPY entrypoint.py /etc/compass/entrypoint.py
ENTRYPOINT ["/etc/compass/entrypoint.sh"]
CMD [""]

