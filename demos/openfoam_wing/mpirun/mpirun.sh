#!/bin/bash

NTASKS="4"
export theImage="library://marcodelapierre/beta/openfoam:v2012"

# this configuration depends on the host
export MPICH_ROOT="/opt/intel/oneapi"
export MPICH_LIBS="$( which mpirun )"
export MPICH_LIBS="${MPICH_LIBS%/bin/mpirun*}/lib"

export SINGULARITY_BINDPATH="$MPICH_ROOT"
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_LIBS:$MPICH_LIBS/libfabric:\$LD_LIBRARY_PATH"

./runAll.sh
