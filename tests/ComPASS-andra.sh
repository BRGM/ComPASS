#!/bin/sh

N_PROC=1

EXE=../build/bin/ComPASS-andra
MESH=meshes/andra_cartesian.msh

OUTPUT=ComPASS-andra
LOG=${OUTPUT}/output.log 

mkdir -p ${OUTPUT}

mpirun -n ${N_PROC} ${EXE} ${MESH} ${LOG} ${OUTPUT}
