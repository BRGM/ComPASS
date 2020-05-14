FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest

# The following two lines are necessary for the configuration of tzdata
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install --yes --no-install-recommends \
    vim  \
    cmake-curses-gui \
    python3-scipy \
    ipython3 \
 && apt-get clean \
 && pip3 install sortedcontainers

# petsc4py installation
# docker seems to be unable to delete /wheel directory, so we just delete its content
# might be linked to https://stackoverflow.com/questions/46824425/docker-is-unable-to-delete-a-directory-when-building-images
RUN pip3 install /wheel/petsc4py-*.whl \
 && rm -rvf /wheel/*.whl

WORKDIR /localfs

COPY entrypoint.sh /etc/compass/entrypoint.sh
COPY entrypoint.py /etc/compass/entrypoint.py
ENTRYPOINT ["/etc/compass/entrypoint.sh"]
CMD [""]
