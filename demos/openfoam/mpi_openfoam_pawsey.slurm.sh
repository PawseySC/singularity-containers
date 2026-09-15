#!/bin/bash --login

#SBATCH --job-name=mpi-openfoam-training
#SBATCH --partition=work
#SBATCH --reservation=ContainersTraining
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00

#--- Load the singularity module (Pawsey's mpi-settings flavour):
module load singularity/4.1.0-mpi

#--- Using user's own image:
export SINGULARITY_IMAGE="$MYSOFTWARE/singularity/images/openfoam--v2012.sif"  #Adapt path and name to the correct ones
echo "Using openfoam singularity image:"
echo "SINGULARITY_IMAGE=$SINGULARITY_IMAGE"

#--- Specific settings for the cluster you are on
#(Check the specific guide of the cluster for additional settings)

#--- Execute pre-processing tools: (Note that in this example, we are using uncollated writing)
#(These pre-processing tools are serial by design)
singularity exec $SINGULARITY_IMAGE blockMesh | tee log.blockMesh
singularity exec $SINGULARITY_IMAGE topoSet | tee log.topoSet
singularity exec $SINGULARITY_IMAGE decomposePar -fileHandler uncollated -force | tee log.topoSet

#--- Execute the solver: (Note that in this example, we are using uncollated writing)
#(Solvers use MPI parallelism by design)
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c 1 \
     singularity exec $SINGULARITY_IMAGE simpleFoam -parallel \
     -fileHandler uncollated | tee log.simpleFoam

#--- Execute post-processing tools: (Note that in this example, we are using uncollated writing)
#(These post-processing tools are serial by design)
singularity exec $SINGULARITY_IMAGE \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar
