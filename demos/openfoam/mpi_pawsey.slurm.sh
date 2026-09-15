#!/bin/bash -l

#SBATCH --job-name=mpi-openfoam
#SBATCH --nodes=1
#SBATCH --reservation=ContainersTraining
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:20:00

SINGULARITY_IMAGE="library://marcodelapierre/beta/openfoam:v2012"

# this configuration depends on the host
module load singularity/4.1.0-mpi
module list


# pre-processing
echo "Running OpenFoam pre-processing steps (these are serial applications)"
srun -n 1 \
  singularity exec $SINGULARITY_IMAGE \
  blockMesh | tee log.blockMesh

srun -n 1 \
  singularity exec $SINGULARITY_IMAGE \
  topoSet | tee log.topoSet

srun -n 1 \
  singularity exec $SINGULARITY_IMAGE \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
srun -n $SLURM_NTASKS \
  singularity exec $SINGULARITY_IMAGE \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
srun -n 1 \
  singularity exec $SINGULARITY_IMAGE \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar
