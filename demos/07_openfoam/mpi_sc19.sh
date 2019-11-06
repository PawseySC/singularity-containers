#!/bin/bash -l


#SBATCH --job-name=mpi
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:20:00
#SBATCH --export=NONE


# this configuration depends on the host
export SINGULARITY_BINDPATH="/opt/mpich/mpich-3.1.4/apps"
export SINGULARITYENV_LD_LIBRARY_PATH="/opt/mpich/mpich-3.1.4/apps/lib"


# pre-processing
srun --export=all -n 1 \
  singularity exec openfoam_v1812.sif \
  blockMesh | tee log.blockMesh

srun --export=all -n 1 \
  singularity exec openfoam_v1812.sif \
  topoSet | tee log.topoSet

srun --export=all -n 1 \
  singularity exec openfoam_v1812.sif \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
srun --export=all -n 2 \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
srun --export=all -n 1 \
  singularity exec openfoam_v1812.sif \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar

