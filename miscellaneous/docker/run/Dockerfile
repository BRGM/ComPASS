From debian:buster-slim
ENV PYTHONPATH=/build/meshtoolsModule:/build/compassModule
RUN apt-get update && apt-get install --yes --no-install-recommends libpython3-dev python3-mpi4py python3-numpy metis libcgal-dev libpetsc-real3.8 && apt-get clean 
WORKDIR /build/
ENTRYPOINT ["/bin/bash"]

