---
title: "Molecular dynamics with GPU containers"
teaching: 5
exercises: 10
questions:
objectives:
- Get started with Nvidia GPU containers for HPC applications
keypoints:
- You can run containerised GPU applications using the flag `--nv`
- By using this flag, Singularity will look for GPU drivers in the host for you, and bind mount them in the container
---


> ## Note
>
> To run exercises from this episode on your own, you'll need a machine with a GPU card and GPU drivers installed.  
> There are examples for both using and not using the Slurm scheduler.
{: .callout}


### Run a molecular dynamics simulation on a GPU with containers

For our example we are going to use Gromacs, a quite popular molecular dynamics package, among the ones that have been optimised to run on GPUs through Nvidia containers.

First, let us cd into `demos/gromacs` and download the container image `nvcr.io/hpc/gromacs:2018.2`:

```
$ cd $TUTO/demos/gromacs
$ singularity pull docker://nvcr.io/hpc/gromacs:2018.2
```
{: .bash}

The current directory has got sample input files picked from the collection of [Gromacs benchmark examples](ftp://ftp.gromacs.org/pub/benchmarks/water_GMX50_bare.tar.gz).  In particular, we're going to use the subset `water-cut1.0_GMX50_bare/1536/`. First let's `gunzip` one of the required input files:

```
$ gunzip conf.gro.gz
```
{: .bash}

Now, from a Singularity perspective, all we need to do to run a GPU application from a container is to add the runtime flag `--nv`.  This will make Singularity look for the Nvidia drivers in the host, and mount them inside the container.

On the host system side, when running GPU applications through Singularity the only requirement consists of the Nvidia driver for the relevant GPU card (the corresponding file is typically called `libcuda.so.<VERSION>` and is located in some library subdirectory of `/usr`).

Do not execute the next two commands, let us just have a look at them.

* Preliminary step
  ```
  $ singularity exec --nv gromacs_2018.2.sif gmx grompp -f pme.mdp
  ```
  {: .bash}

* Production step
  ```
  $ singularity exec --nv gromacs_2018.2.sif gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
  ```
  {: .bash}

GPU resources are usually made available in HPC systems through schedulers, to which Singularity natively and transparently interfaces.  So, for instance let us have a look in the current directory at the Slurm batch script called `gpu.sh`, which is suitable for running at Pawsey on *Zeus* or *Topaz*:

```
#!/bin/bash -l

#SBATCH --job-name=gpu
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00

module load singularity

# run Gromacs preliminary step with container
srun singularity exec --nv gromacs_2018.2.sif \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun singularity exec --nv gromacs_2018.2.sif \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
```
{: .bash}

Basically, we have just combined the Slurm command `srun` with `singularity exec <..>` (similar to what we did in the episode on MPI).  We can submit the script with:

```
$ sbatch gpu.sh
```
{: .bash}


> ## Nvidia GPU Cloud, a useful resource
>
>The GPU manufacturer *Nvidia* has a dedicated web registry for container images, shipping GPU optimised applications: <https://ngc.nvidia.com>.
>
>To browse this registry, you'll need a free account.  Go to <https://ngc.nvidia.com>, complete the procedure to **Create an account**, then **Sign in** (currently both >options are available on the top right corner of the page).
>
>You can browse the available containerised packages through the various category boxes.  *E.g.* click on the **High Performance Computing** box, then click on the >**Gromacs** one.
{: .callout}
