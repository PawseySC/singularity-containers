---
title: "Computational Fluid Dynamics with MPI containers"
teaching: 10
exercises: 10
questions:
objectives:
- Learn the steps required to configure and run MPI applications from a container
keypoints:
- You need to build your application in the container with an MPI version which is ABI compatible with MPI libraries in the host
- Appropriate environment variables and bind mounts are required at runtime to make the most out of MPI applications (sys admins are your best friends)
- Singularity transparently interfaces with HPC schedulers such as Slurm
---


### Let's run OpenFoam in a container!

We're going to start this episode with actually running a practical example. We'll discuss the way this all works later on.  
We're using OpenFoam, a widely popular package for Computational Fluid Dynamics simulations which is able to massively scale in parallel architectures up to thousands of processes, using the MPI library.  
The sample inputs come straight from the OpenFoam installation tree, namely `$FOAM_TUTORIALS/incompressible/pimpleFoam/LES/periodicHill/steadyState/`.

First, let us cd into the demo directory, and download an appropriate container image:

```
$ cd $SC19/demos/07_openfoam
$ singularity pull library://marcodelapierre/beta/openfoam:v1812
```
{: .bash}

Now, let us use the Slurm scheduler to submit the job script `mpi_sc19.sh`, that will run the sample simulation:

```
$ sbatch mpi_sc19.sh
```
{: .bash}

The run will take a couple of minutes. When it's finished, the directory contents will look a bit like this one:

```
$ ls -ltr
```
{: .bash}

```
total 1121656
drwxr-sr-x  2 mdelapierre pawsey0001       4096 Nov  5 15:45 0
-rwxr-x---+ 1 mdelapierre pawsey0001 1148433339 Nov  6 23:54 openfoam_v1812.sif
-rw-rw----+ 1 mdelapierre pawsey0001        927 Nov  7 00:02 update-settings.sh
drwxr-sr-x  2 mdelapierre pawsey0001       4096 Nov  7 00:02 system
-rw-rw----+ 1 mdelapierre pawsey0001        775 Nov  7 00:25 mpi.sh
drwxrws---+ 4 mdelapierre pawsey0001       4096 Nov  7 00:27 dynamicCode
drwxr-sr-x  3 mdelapierre pawsey0001       4096 Nov  7 00:27 constant
-rw-rw----+ 1 mdelapierre pawsey0001       3594 Nov  7 00:27 log.blockMesh
-rw-rw----+ 1 mdelapierre pawsey0001       1948 Nov  7 00:27 log.topoSet
-rw-rw----+ 1 mdelapierre pawsey0001       2311 Nov  7 00:28 log.decomposePar
drwxrws---+ 8 mdelapierre pawsey0001       4096 Nov  7 00:29 processor1
drwxrws---+ 8 mdelapierre pawsey0001       4096 Nov  7 00:29 processor0
-rw-rw----+ 1 mdelapierre pawsey0001      18573 Nov  7 00:29 log.simpleFoam
-rw-rw----+ 1 mdelapierre pawsey0001      28224 Nov  7 00:29 slurm-4198976.out
-rw-rw----+ 1 mdelapierre pawsey0001       1540 Nov  7 00:29 log.reconstructPar
drwxrws---+ 3 mdelapierre pawsey0001       4096 Nov  7 00:29 20
```
{: .output}

We ran using 2 MPI processes, who created outputs in the directories `processor0` and `processor1`, respectively. The final reconstruction creates results in the directory `20` (which stands for the 20th and last simulation step in this very short demo run).

What has just happened?


### A batch script for MPI applications with containers

Let's have a look at the content of the script (`mpi_sc19.sh`) we executed through the scheduler:

```
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
```
{: .bash}

