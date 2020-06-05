#!/bin/bash -l

#SBATCH --job-name=pull_images
#SBATCH --ntasks=1
#SBATCH --time=06:00:00

# executing a dummy command just to cause images to be downloaded in the cache

# day 1
singularity exec docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4 echo ciao

# day 2
singularity exec docker://rocker/tidyverse:3.6.1 echo ciao
singularity exec docker://jupyter/datascience-notebook:latest echo ciao

# day 3
singularity exec library://marcodelapierre/beta/openfoam:v1812 echo ciao
singularity exec docker://nvcr.io/hpc/gromacs:2018.2 echo ciao
singularity exec docker://trinityrnaseq/trinityrnaseq:2.8.6 echo ciao

# bonus
singularity exec docker://nextflow/rnaseq-nf:latest echo ciao
singularity exec docker://marcodelapierre/gnuplot:5.2.2_4 echo ciao
# codimd: nope, no docker on HPC
