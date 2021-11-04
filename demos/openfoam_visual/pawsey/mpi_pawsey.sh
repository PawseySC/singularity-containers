#!/bin/bash -l

#SBATCH --job-name=mpi
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:20:00

unset SBATCH_EXPORT
export theImage="library://marcodelapierre/beta/openfoam:v2012"

# this configuration depends on the host
module unload xalt
module load singularity

./runAll.sh
