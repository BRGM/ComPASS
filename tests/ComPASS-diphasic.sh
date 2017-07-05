#!/bin/sh

N_PROC=1

EXE=../build/bin/ComPASS-diphasic
MESH=meshes/gallery.msh

OUTPUT=ComPASS-diphasic
LOG=${OUTPUT}/output.log 

mkdir -p ${OUTPUT}

mpirun -n ${N_PROC} ${EXE} ${MESH} ${LOG} ${OUTPUT}
