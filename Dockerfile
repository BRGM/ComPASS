FROM registry.gitlab.inria.fr/charms/compass/build-environment:latest AS build

RUN mkdir -p /source
COPY ./ /source/ComPASS
WORKDIR /build

# If you want to compile only specific physics comment the following line and
# use the next run command
RUN CC=mpicc cmake ../source/ComPASS && make -j `nproc`
#RUN CC=mpicc cmake ../source/ComPASS \
#    -DComPASS_WITH_diphasic_PHYSICS=0, \
#    -DComPASS_WITH_linear_water_PHYSICS=1, \
#    -DComPASS_WITH_linear_liquid_water_PHYSICS=0, \
#    -DComPASS_WITH_liquid_water_PHYSICS=0, \
#    -DComPASS_WITH_water2ph_PHYSICS=0, \
#    -DComPASS_WITH_water_with_tracer_PHYSICS=0 \
# && make -j `nproc`


FROM registry.gitlab.inria.fr/charms/compass/run-environment:latest

RUN mkdir -p /build/ComPASS
COPY --from=build /source/ComPASS/python /build/ComPASS/python
ENV PYTHONPATH=/build/ComPASS/python

WORKDIR /localfs
