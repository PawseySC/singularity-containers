#!/bin/bash -l

#SBATCH --job-name=mpi
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:20:00

image="library://marcodelapierre/beta/openfoam:v1812"

# this configuration depends on the host
module unload xalt
module load singularity


# pre-processing
srun -n 1 \
  singularity exec $image \
  blockMesh | tee log.blockMesh

srun -n 1 \
  singularity exec $image \
  topoSet | tee log.topoSet

srun -n 1 \
  singularity exec $image \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
srun -n $SLURM_NTASKS \
  singularity exec $image \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
srun -n 1 \
  singularity exec $image \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar

