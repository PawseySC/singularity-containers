#!/bin/bash

NTASKS="4"
export theImage="library://marcodelapierre/beta/openfoam:v2012"

# this configuration depends on the host
export MPICH_ROOT="/opt/mpich/mpich-3.1.4/apps"

export SINGULARITY_BINDPATH="$MPICH_ROOT"
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_ROOT/lib:\$LD_LIBRARY_PATH"

./runAll.sh
