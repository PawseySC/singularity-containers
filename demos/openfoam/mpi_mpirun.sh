#!/bin/bash

NTASKS="2"

# this configuration depends on the host
export SINGULARITY_BINDPATH="/opt/mpich/mpich-3.1.4/apps"
export SINGULARITYENV_LD_LIBRARY_PATH="/opt/mpich/mpich-3.1.4/apps/lib"


# pre-processing
singularity exec openfoam_v1812.sif \
  blockMesh | tee log.blockMesh

singularity exec openfoam_v1812.sif \
  topoSet | tee log.topoSet

singularity exec openfoam_v1812.sif \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
mpirun -n $NTASKS \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
singularity exec openfoam_v1812.sif \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar

