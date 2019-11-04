---
layout: page
title: Setup
root: .
---

The main requirement for this workshop is a personal computer with a web browser and a command line shell program.  
These will allow you to follow the online materials and to login to a facility with the required software stack.

Options: 

* **SC19 Workshop**: you will will need username and password to access the training VMs, which will be provided by instructors.

* **Zeus @Pawsey**: if you have access, Singularity can be loaded with `module load singularity/3.3.0`. MPI libraries are configured properly, and GPU applications can be run on the Slurm partition `gpuq`.

* **Nimbus Cloud @Pawsey**: if you have access, both Singularity and Docker are preinstalled in the *Ubuntu Pawsey* base image.

* **BYO device**: if you have a Linux box, you can install the required software yourself (might take a while):

  * Essential (core of the tutorial)
    - Singularity
  
  * Desirable (to run all Singularity examples)
    - Slurm scheduler
    - MPICH library
    - Nvidia GPU driver (GPU card required)
  
  * Optional (extra applications for last two episodes)
    - Docker
    - HPCCM
    - Podman
    - Sarus
