#!/bin/bash -l

#SBATCH --job-name=gpu
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00

image="docker://nvcr.io/hpc/gromacs:2018.2"
module load singularity

# run Gromacs preliminary step with container
srun singularity exec --nv $image \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun singularity exec --nv $image \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
