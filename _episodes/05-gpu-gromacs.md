---
title: "Molecular dynamics with GPU containers"
teaching: 5
exercises: 10
questions:
objectives:
- Get started with Nvidia GPU containers for HPC applications
keypoints:
- You can run containerised GPU applications using the flag `--nv`
- Singularity transparently interfaces with HPC schedulers such as Slurm
---


### Nvidia GPU Cloud

The GPU manufacturer *Nvidia* has a dedicated web registry for container images, shipping GPU optimised applications: <https://ngc.nvidia.com>.

To browse this registry, you'll need a free account. Go to <https://ngc.nvidia.com>, complete the procedure to **Create an account**, then **Sign in** (currently both options are available on the top right corner of the page).

You can browse the available containerised packages through the various category boxes. E.g. click on the **High Performance Computing** box, then click on the **Gromacs** one. The page will briefly discuss the code, with instructions on how to pull and run the container.


### Run a molecular dynamics simulation on a GPU with containers

For our example we are going to use Gromacs, a quite popular molecular dynamics package, among the ones that have been optimised to run on GPUs through Nvidia containers.

First, let us cd into `demos/05_gromacs`:

```
$ cd $SC19/demos/05_gromacs
```
{: .bash}

> Now, try and `pull` the container `nvcr.io/hpc/gromacs:2018.2`. **Hint**: container images in this registry are in Docker format, so you will need to prepend the prefix `docker://`.
> 
> > ## Solution
> > 
> > ```
> > $ singularity pull docker://nvcr.io/hpc/gromacs:2018.2'
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

The current directory has got some sample input files picked from the collection of [Gromacs benchmark examples](ftp://ftp.gromacs.org/pub/benchmarks/water_GMX50_bare.tar.gz).


THIS NEEDS TO BE PREPARED
```
$ wget ftp://ftp.gromacs.org/pub/benchmarks/water_GMX50_bare.tar.gz
$ tar xzf water_GMX50_bare.tar.gz
$ cp water-cut1.0_GMX50_bare/1536/* .
```
{: .bash}

Now, from a Singularity perspective, all we need to do to run a GPU application from a container is to add the runtime flag `--nv`. This will make Singularity look for the Nvidia drivers in the host, and mount them inside the container.

GPU resources are usually made available in HPC systems through schedulers, to which Singularity natively and transparently interfaces. This is how a batch script ..


Similar to the GPU machine learning example in a previous episode, only minor modifications are required in the SLURM script to run on GPUs:

* the GPU partition on Zeus is called `gpuq`, this is in substitution for the `workq` we've used in CPU jobs;
* we need to set an additional SBATCH flag, `--gres=gpu:1`, to request use of a GPU.

Use your favourite text editor to create a script `gpu.sh`:

```
#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --job-name=gpu

module load shifter

# run Gromacs preliminary step with container
srun --export=all shifter run nvcr.io/hpc/gromacs:2018.2 \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun --export=all shifter run nvcr.io/hpc/gromacs:2018.2 \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
```
{: .bash}

The script is ready for submission:

```
$ sbatch --reservation <your-pawsey-reservation> gpu.sh
```
{: .bash}

