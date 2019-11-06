#!/bin/bash -l

#SBATCH --job-name=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --export=NONE

# run Gromacs preliminary step with container
srun --export=all \
    singularity exec --nv docker://nvcr.io/hpc/gromacs:2018.2 \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun --export=all \
    singularity exec --nv docker://nvcr.io/hpc/gromacs:2018.2 \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
