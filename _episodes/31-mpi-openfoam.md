---
title: "Computational Fluid Dynamics with MPI containers"
teaching: 10
exercises: 10
questions:
objectives:
- Discuss the steps required to configure and run MPI applications from a container
- Discuss the performance of parallel applications inside containers *versus* regular runs
keypoints:
- You need to build your application in the container with an MPI version which is ABI compatible with MPI libraries in the host
- Appropriate environment variables and bind mounts are required at runtime to make the most out of MPI applications (sys admins can help)
- Singularity interfaces almost transparently with HPC schedulers such as Slurm
- MPI performance of containerised applications almost coincide with those of a native run
---


### Let's run OpenFoam in a container!

We're going to start this episode with actually running a practical example, and then discuss the way this all works later on.  
We're using OpenFoam, a widely popular package for Computational Fluid Dynamics simulations, which is able to massively scale in parallel architectures up to thousands of processes, by leveraging an MPI library.  
The sample inputs come straight from the OpenFoam installation tree, namely `$FOAM_TUTORIALS/incompressible/pimpleFoam/LES/periodicHill/steadyState/`.

First, let us cd into the demo directory and download the OpenFoam container image:

```
$ cd $TUTO/demos/openfoam
$ singularity pull library://marcodelapierre/beta/openfoam:v1812
```
{: .bash}

Now, let us run the sample simulation with:

```
$ ./mpi_mpirun.sh
```
{: .bash}

**In alternative**, if you're running this example on Pawsey systems (*e.g.* Magnus or Zeus), achieve the same result by using the Slurm scheduler to submit the job script `mpi_pawsey.sh`:

```
$ sbatch mpi_pawsey.sh
```
{: .bash}

The run will take a couple of minutes. When it's finished, the directory contents will look a bit like this one:

```
$ ls -ltr
```
{: .bash}

```
total 80
-rwxr-xr-x 1 user000 tutorial  1304 Nov 16 17:36 update-settings.sh
drwxr-xr-x 2 user000 tutorial   141 Nov 16 17:36 system
-rw-r--r-- 1 user000 tutorial   871 Nov 16 17:36 mpi_pawsey.sh
-rwxr-xr-x 1 user000 tutorial   789 Nov 16 17:36 mpi_mpirun.sh
drwxr-xr-x 2 user000 tutorial    59 Nov 16 17:36 0
drwxr-xr-x 4 user000 tutorial    72 Nov 16 22:45 dynamicCode
drwxr-xr-x 3 user000 tutorial    77 Nov 16 22:45 constant
-rw-rw-r-- 1 user000 tutorial  3493 Nov 16 22:45 log.blockMesh
-rw-rw-r-- 1 user000 tutorial  1937 Nov 16 22:45 log.topoSet
-rw-rw-r-- 1 user000 tutorial  2300 Nov 16 22:45 log.decomposePar
drwxr-xr-x 8 user000 tutorial    70 Nov 16 22:47 processor1
drwxr-xr-x 8 user000 tutorial    70 Nov 16 22:47 processor0
-rw-rw-r-- 1 user000 tutorial 18569 Nov 16 22:47 log.simpleFoam
drwxr-xr-x 3 user000 tutorial    76 Nov 16 22:47 20
-rw-r--r-- 1 user000 tutorial 28617 Nov 16 22:47 slurm-10.out
-rw-rw-r-- 1 user000 tutorial  1529 Nov 16 22:47 log.reconstructPar
```
{: .output}

We ran using *2 MPI* processes, who created outputs in the directories `processor0` and `processor1`, respectively.  The final reconstruction creates results in the directory `20` (which stands for the *20th* and last simulation step in this very short demo run).

What has just happened?


### A batch script for MPI applications with containers

Let's have a look at the content of the script `mpi_mpirun.sh`:

```
#!/bin/bash

NTASKS="2"

# this configuration depends on the host
export MPICH_ROOT="/opt/mpich/mpich-3.1.4/apps"

export SINGULARITY_BINDPATH="$MPICH_ROOT"
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_ROOT/lib"


# pre-processing
singularity exec openfoam_v1812.sif \
  blockMesh | tee log.blockMesh

singularity exec openfoam_v1812.sif \
  topoSet | tee log.topoSet

singularity exec openfoam_v1812.sif \
  decomposePar -fileHandler uncollated | tee log.decomposePar


# run OpenFoam with MPI
mpirun -n $NTASKS \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam


# post-processing
singularity exec openfoam_v1812.sif \
  reconstructPar -latestTime -fileHandler uncollated | tee log.reconstructPar
```
{: .bash}


### How does Singularity interplay with the MPI launcher?

We'll comment on the environment variable definitions soon, now let's focus on the set of commands that make the simulation happen.

In particular, the fourth command is the only one using multiple processors through MPI:

```
mpirun -n $NTASKS \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam
```
{: .bash}

Here, `mpirun` is the MPI launcher, *i.e.* the tool that is in charge for spawning the multiple MPI processes that will make the workflow run in parallel.  
Note how `singularity` can be executed through the launcher as any other application would.

Under the hood, the MPI processes outside of the container (spawned by `mpirun`) will work in tandem with the containerized MPI code to instantiate the job.  
There are a few implications here...


### Requirements for the MPI + container combo

Let's discuss what the above mentioned implications are.

* A host MPI installation must be present to spawn the MPI processes.

* An MPI installation is required in the container, to compile the application.  Also, during build the application must be linked *dynamically* to the MPI libraries, so as to have the capability of using the host ones at runtime.  Note how dynamic linking is typically the default behaviour on Linux systems.  
A specific section of the recipe file needs to take care of this, or in alternative the base image for the recipe needs to have the MPI libraries.  Either way, if we take the example of a *def file* for the *MPICH* flavour of MPI, the code would look like:

```
%post

[..]

MPICH_VERSION="3.1.4"
MPICH_CONFIGURE_OPTIONS="--enable-fast=all,O3 --prefix=/usr"

mkdir -p /tmp/mpich-build
cd /tmp/mpich-build

wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz
tar xvzf mpich-${MPICH_VERSION}.tar.gz

cd mpich-${MPICH_VERSION}

./configure ${MPICH_CONFIGURE_OPTIONS}
make
make install

ldconfig

[..]
```
{: .bash}


> ## Base MPI image at Pawsey
>
> Pawsey maintains an MPICH base image at [pawsey/mpich-base](https://hub.docker.com/r/pawsey/mpich-base).  
> At the moment, only a Docker image is provided, which of course can also be used by Singularity.
{: .callout}


* The container and host MPI installations need to be *ABI* (Application Binary Interface) *compatible*. This is because the application in the container is built with the former but runs with the latter.  
At present, there are just two families of MPI implementations, not ABI compatible with each other: MPICH (with IntelMPI and MVAPICH) and OpenMPI.  
If you anticipate your application will run in systems with non ABI compatible libraries, you will need to build variants of the image for the two MPI families.


> ## MPI implementations at Pawsey
>
> At present, all Pawsey systems have installed at least one MPICH ABI compatible implementation: CrayMPICH on the Crays (*Magnus* and *Galaxy), IntelMPI on *Zeus* and *Topaz*.  Therefore, MPICH is the recommended MPI library to install in container images.  
> Zeus and Topaz also have OpenMPI, so images built over this MPI family can run in these clusters, upon appropriate configuration of the shell environment (see below).
{: .callout}


* Bind mounts and environment variables need to be setup so that the containerised MPI application can use the host MPI libraries at runtime.  Bind mounts can be configured by the administrators, or set up through variables. We're discussing the latter way here.  
In the current example we have:

```
export MPICH_ROOT="/opt/mpich/mpich-3.1.4/apps"

export SINGULARITY_BINDPATH="$MPICH_ROOT"
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_ROOT/lib"
```
{: .bash}

Here, `SINGULARITY_BINDPATH` bind mounts the host path where the MPI installation is (MPICH in this case).  
The second variable, SINGULARITYENV_LD_LIBRARY_PATH, ensures that at runtime the container's `LD_LIBRARY_PATH` has the path to the MPICH libraries.

> ## Interconnect libraries and containers
>
> If the HPC system you're using has high speed interconnect infrastructure, than it will also have some system libraries to handle that at the application level.  These libraries will need to be exposed to the containers, too, similar to the MPI libraries, to ensure maximum performance are achieved.  
> This can be a challenging task for a user, as it requires knowing details on the installed software stack.  System administrators should be able to assist in this regard.
{: .callout}

> ## Singularity environment variables at Pawsey
>
> In all Pawsey systems, the Singularity module sets up all of the required variables for MPI and interconnect libraries.  So this will do the job:
>
> ```
> $ module load singularity
> ```
> {: .bash}
{: .callout}


### Singularity interface to Slurm

Now, if we have a look at the script variant for the Slurm scheduler, `mpi_pawsey.sh`, we'll see the key difference is that every OpenFoam command is executed via `srun`:

```
srun -n $SLURM_NTASKS \
  singularity exec openfoam_v1812.sif \
  simpleFoam -fileHandler uncollated -parallel | tee log.simpleFoam
```
{: .bash}

`srun` is the Slurm wrapper for the MPI launcher, `mpirun`.  Other schedulers will require a different command.  
In practice, all we had to do was to replace `mpirun` with `srun`.  This is because Singularity implements a native interface to schedulers, so it can be executed through `srun` as other packages would.

Note in the script how, when using schedulers, it is good practice to execute all application commands through `srun`, even those that only use one core.  


### MPI performance: container *vs* bare metal

What's the performance overhead in running an MPI application through containers?

Well, the benchmark figures just below reveal it's quite small...good news!

<!-- ![OSU bandwidth test]({{ page.root }}/fig/OSU_Bandwidth.png) -->
<img src="{{ page.root }}/fig/OSU_Bandwidth.png" alt="OSU bandwidth test" width="651" height="489"/>

<!-- ![OSU point-to-point latency test]({{ page.root }}/fig/OSU_Latency_P2P.png) -->
<img src="{{ page.root }}/fig/OSU_Latency_P2P.png" alt="OSU point-to-point latency test" width="651" height="489"/>

<!-- ![OSU collective latency test]({{ page.root }}/fig/OSU_Latency_Coll.png) -->
<img src="{{ page.root }}/fig/OSU_Latency_Coll.png" alt="OSU collective latency test" width="651" height="489"/>
