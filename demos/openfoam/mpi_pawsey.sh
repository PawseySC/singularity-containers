#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --job-name=mpi
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:20:00


# this configuration depends on the host
module load singularity



# pre-processing
srun -n 1 \
  singularity exec openfoam_v1812.sif \
  blockMesh | tee log.blockMesh

srun -n 1 \
  singularity exec openfoam_v1812.sif \
  topoSet | tee log.topoSet

srun -n 1 \
  singularity exec openfoam_v1812.sif \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
srun -n $SLURM_NTASKS \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
srun -n 1 \
  singularity exec openfoam_v1812.sif \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar

