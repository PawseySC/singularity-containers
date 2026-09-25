#!/bin/bash --login

#SBATCH --job-name=mpi-openfoam-training
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

#--- Load the singularity module (Pawsey's mpi-settings flavour):
module load singularity/4.1.0-mpi

#--- Using user's own image:
export SINGULARITY_IMAGE="$MYSOFTWARE/singularity/images/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif"  #Adapt path and name to the correct ones
echo "Using openfoam singularity image:"
echo "SINGULARITY_IMAGE=$SINGULARITY_IMAGE"

#--- Prepare the case directory:
caseDir=periodicPlaneChannel
cd $caseDir
rm -rf 0
cp -r 0.orig 0

#--- Specific settings for the cluster you are on
#(Check the specific guide of the cluster for additional settings)

#--- Automating the list of IORANKS for collated fileHandler
echo "Setting the grouping ratio for collated fileHandling"
nProcs=$SLURM_NTASKS #Number of total processors in decomposition for this case
mGroup=4             #Size of the groups for collated fileHandling (32 is the initial recommendation for Setonix)
of_ioRanks="0"
iC=$mGroup
while [ $iC -lt $nProcs ]; do
   of_ioRanks="$of_ioRanks $iC"
   ((iC += $mGroup))
done
export FOAM_IORANKS="("${of_ioRanks}")"
echo "FOAM_IORANKS=$FOAM_IORANKS"

#--- Execute pre-processing tools:
#(These pre-processing tools are serial by design)
singularity exec $SINGULARITY_IMAGE blockMesh | tee log.blockMesh
singularity exec $SINGULARITY_IMAGE renumberMesh -overwrite -constant | tee log.renumberMesh
singularity exec $SINGULARITY_IMAGE decomposePar -cellDist -force | tee log.decomposePar

#--- Execute the parallel solver:
#(Solvers use MPI parallelism by design)
srun -N $SLURM_JOB_NUM_NODES -n $nProcs -c 1 \
  singularity exec $SINGULARITY_IMAGE pimpleFoam -parallel | tee log.pimpleFoam

#--- Execute post-processing tools:
#(These post-processing tools are serial by design)
singularity exec $SINGULARITY_IMAGE reconstructPar -latestTime | tee log.reconstructPar
singularity exec $SINGULARITY_IMAGE postChannel -latestTime | tee log.postChannel

#--- Final commands
echo "OpenFOAM script has reached the end"
