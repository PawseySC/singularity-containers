#!/bin/bash

NTASKS="2"
image="library://marcodelapierre/beta/openfoam:v1812"

# this configuration depends on the host
export MPICH_ROOT="/opt/mpich/mpich-3.1.4/apps"

export SINGULARITY_BINDPATH="$MPICH_ROOT"
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_ROOT/lib"



# pre-processing
singularity exec $image \
  blockMesh | tee log.blockMesh

singularity exec $image \
  topoSet | tee log.topoSet

singularity exec $image \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
mpirun -n $NTASKS \
  singularity exec $image \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
singularity exec $image \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar

