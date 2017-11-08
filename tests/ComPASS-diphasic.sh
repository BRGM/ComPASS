#!/bin/sh

N_PROC=1

EXE=../build/bin/ComPASS-diphasic
MESH=meshes/gallery.msh

OUTPUT=ComPASS-diphasic
LOG=${OUTPUT}/output.log 

rm -r proc.*
rm -r ${OUTPUT}/{time_*,*data.pvd}
mkdir -p ${OUTPUT}

mpirun -n ${N_PROC} ${EXE} ${MESH} ${LOG} ${OUTPUT}
