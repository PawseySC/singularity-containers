#!/bin/bash -l

#SBATCH --job-name=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00

# run Gromacs preliminary step with container
srun singularity exec --nv $SIFPATH/gromacs_2018.2.sif \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun singularity exec --nv $SIFPATH/gromacs_2018.2.sif \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
