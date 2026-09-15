---
title: "Using MPI-capable containers (with OpenFOAM as an example)"
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

### Are you running on a shared HPC system?

If you're running this tutorial on a shared system (*e.g.* Setonix at Pawsey), you should use one of the compute nodes rather than the login node.  You can set this up by using an interactive scheduler allocation, for instance on Setonix with Slurm (do this if you are not in an `salloc` interactive session yet):

```
$ salloc -N 1 -n 1 -c 8 --reservation=ContainersTraining -t 4:00:00
```
{: .source}

```text
salloc: Granted job allocation 3453895
salloc: Waiting for resource configuration
salloc: Nodes nid000152 are ready for job
```
{: .output}

### Get ready for the hands-on

Before we start, let us ensure we have the required files to run the tutorials.

If you haven't done so already, move to a suitable working directory and download the following GitHub repository. On Pawsey systems, use your scratch directory; on other HPC or cloud systems, use the equivalent working directory recommended by the system administrators.

```bash
$ cd "$MYSCRATCH"    # On Pawsey systems
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
$ cd "$TUTO"
```
{: .source}

Now `cd` to the working directory. In this case:
```bash
$ cd demos/openfoam
$ pwd
```

The working directory should be something like:
```text
/path/to/scratch/singularity-containers/demos/openfoam
```
{: .output}

Load the singularity module (in this case, Pawsey's mpi-ready flavour):

```bash
$ module load singularity/4.1.0-mpi
```
{: .source}

### Choose an OpenFOAM image provided by Pawsey in the quay.io registry

In your own webbrowser, go to `https://quay.io/pawsey`.

Pawsey provides several images useful for maby different research areas. From there you should be able to see a repository named `pawsey/openfoam`. Click on it and then in the tags icon (second icon top to bottom on the left side of the screen). (**Do not confuse with the other flavour named `openfoam-org`.**)

You should be able to see the openfoam image with the tag: `v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04`. You should be able click the fetch tag ico on the right and then choose the format `Docker Pull (by tag)`. Then, copy just the tag (**NOT THE DOCKER COMMAND**) and close the prompt.

### Pull the image to be used into your personal library

In the interactive shell terminal connected to the salloc session in Setonix, define your personal library directory (and create the directory if not done yet):

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```

Pull the image to use for this example into your personal library directory:

```bash
$ singularity pull \
  "${MY_LOCAL_LIBRARY}/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif" \
   docker://quay.io/pawsey/openfoam:v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04
```
{: .source}

This may take some time, so you can go for a coffee in the meantime.

### Prepare the tutorial starting with a copy from the image into the working directory in the host

OpenFOAM counts with several examples (tutorials) that can be used as starting point for your research or learning process. Selection of a tutorial is usually performed in an interactive shell inside the container. First, start the interactive shell (you will notice that prompt will change to `Singularity>`):

```bash
$ SINGULARITY_IMAGE="${MY_LOCAL_LIBRARY}/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif"
$ singularity shell "$SINGULARITY_IMAGE"
```
{: .source}

```text
Singularity>
```
{: .output}

Check that, at entrance, the working directory is the same from which the singularity shell was invoked (and save that path in a variable):

```bash
Singularity> pwd
Singularity> HOST_WORKING_DIR="$(pwd)"
Singularity> echo "$HOST_WORKING_DIR"
```
{: .source}

OpenFOAM defines several environment variables to make your life easier. In this case you can explore the available tutorials under the directory defined by `FOAM_TUTORIALS`. Go into that directory using `cd` and list its content:

```bash
Singularity> echo "$FOAM_TUTORIALS"
Singularity> cd "$FOAM_TUTORIALS"
Singularity> pwd
Singularity> ls
```
{: .source}

```text
Allclean    DNS		compressible	  finiteArea	  mesh		 resources
Allcollect  IO		discreteMethods   heatTransfer	  modules	 stressAnalysis
Allrun	    basic	electromagnetics  incompressible  multiphase	 verificationAndValidation
Alltest     combustion	financial	  lagrangian	  preProcessing
```
{: .output}

(At this point we will assume that the selected tutorial is `$FOAM_TUTORIALS/incompressible/pimpleFoam/LES/periodicPlaneChannel`.)

Once the tutorial have been selected, copy the case directory into the working directory in the host (we previously saved that path in `HOST_WORKING_DIR` variable):

```bash
Singularity> cp -r "${FOAM_TUTORIALS}/incompressible/pimpleFoam/LES/periodicPlaneChannel" "$HOST_WORKING_DIR"
```
{: .source}

> ## Alternative: copy the tutorial without opening an interactive shell
>
> The same steps can be performed directly from the host by using `singularity exec`. First, inspect the location and contents of the OpenFOAM tutorials directory:
>
> ```bash
> $ singularity exec "$SINGULARITY_IMAGE" \
>     bash -c 'find $FOAM_TUTORIALS -iname "*PlaneChannel*"'
> ```
> {: .source}
>
> The resulting list would look something like this:
>
> ```text
> /opt/OpenFOAM/OpenFOAM-v2606/tutorials/incompressible/pimpleFoam/LES/periodicPlaneChannel
> /opt/OpenFOAM/OpenFOAM-v2606/tutorials/incompressible/pimpleFoam/LES/planeChannel
> /opt/OpenFOAM/OpenFOAM-v2606/tutorials/verificationAndValidation/turbulenceModels/planeChannel
> /opt/OpenFOAM/OpenFOAM-v2606/tutorials/verificationAndValidation/turbulentInflow/oneCellThickPlaneChannel
> ```
> {: .output}
>
> Then copy the selected tutorial into the current working directory on the host:
>
> ```bash
> $ singularity exec "$SINGULARITY_IMAGE" \
>     bash -c 'cp -r "$FOAM_TUTORIALS/incompressible/pimpleFoam/LES/periodicPlaneChannel" "$PWD"'
> ```
> {: .source}
>
> The command is passed through `bash -c` so that `$FOAM_TUTORIALS` and `$PWD` are expanded inside the container. Singularity makes the host working directory available inside the container at the same path, so the copied `periodicPlaneChannel` directory appears in the directory from which the command was run.
{: .solution}

Now update the default settings of the OpenFOAM tutorial to preferred settings for this tutorial. This is done by the `update-settings.sh` script.

```bash
$ ./update-settings.sh
```
{: .source}

> ## If curious: inspect the changes
>
> If you are curious about what was updated, check the differences of the OpenFOAM dictionaries in `periodicPlaneChannel/system` subdirectory with respect to their original settings now copied to `*.original.00` files.
{: .solution}


### Run the MPI containerised application in an HPC cluster!

Submit the Slurm job script:

```bash
$ sbatch --reservation=ContainersTraining mpi_openfoam_pawsey.slurm.sh
```
{: .source}

Check the execution status of the job with:

```bash
$ squeue --me
```
{: .source}

```text
JOBID        USER ACCOUNT             NAME EXEC_HOST ST  REASON START_TIME   END_TIME  TIME_LEFT NODES   PRIORITY     QOS
48927321 course01 courses   mpi-openfoam-t nid002604  R    None 18:40:12     19:00:12      19:45     1      75246  normal
```
{: .output}

Output of the job could be monitored in the slurm output file. For example `slurm-48927321.out`. Use `tail -f` to have a live update of the progress (the name of your file will be different):

```bash
$ tail -f slurm-48927321.out
```
{: .source}

Exit the display with `<Ctrl>-C`

Once the job has finished, the results in the case directory will look something like this:

```bash
$ ls -ltr periodicPlaneChannel
```
{: .source}

```
-rwxr-xr-x   1 courses01 courses     915 Sep 15 19:05 Allrun
-rwxr-xr-x   1 courses01 courses     340 Sep 15 19:05 Allclean
drwxr-sr-x   2 courses01 courses    4096 Sep 15 19:05 0.orig
drwxr-sr-x   2 courses01 courses    4096 Sep 15 19:06 system
-rw-r--r--   1 courses01 courses    3299 Sep 15 19:08 log.blockMesh
-rw-r--r--   1 courses01 courses    2167 Sep 15 19:08 log.renumberMesh
drwxr-sr-x   3 courses01 courses    4096 Sep 15 19:08 constant
drwxr-sr-x   2 courses01 courses    4096 Sep 15 19:08 0
-rw-r--r--   1 courses01 courses    5782 Sep 15 19:08 log.decomposePar
drwxr-sr-x 104 courses01 courses    4096 Sep 15 19:08 processors8_4-7
drwxr-sr-x 104 courses01 courses    4096 Sep 15 19:08 processors8_0-3
-rw-r--r--   1 courses01 courses 1160982 Sep 15 19:08 log.pimpleFoam
-rw-r--r--   1 courses01 courses    2092 Sep 15 19:09 log.reconstructPar
drwxr-sr-x   3 courses01 courses    4096 Sep 15 19:09 200
-rw-r--r--   1 courses01 courses    1763 Sep 15 19:09 log.postChannel
drwxr-sr-x   3 courses01 courses    4096 Sep 15 19:13 graphs
```
{: .output}

We ran using *8 MPI* processes, who created outputs in the directories `processors8_0-3` and `processors8_4-7`.  The final reconstruction has only been applied to the latest time and creates results in the directory `200` (which stands for simulation time `200`).

### A batch script for MPI applications with containers

Let's get back to the directory path for the first example:

```
$ cd $TUTO/demos/openfoam
```
{: .bash}

and have a look at the content of the script `mpi_openfoam_pawsey.slurm.sh`:

```bash
#!/bin/bash --login

#SBATCH --job-name=mpi-openfoam-training
#SBATCH --partition=work
#SBATCH --reservation=ContainersTraining
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
while [ $iC -le $nProcs ]; do
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
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c 1 \
  singularity exec $SINGULARITY_IMAGE pimpleFoam -parallel | tee log.pimpleFoam

#--- Execute post-processing tools:
#(These post-processing tools are serial by design)
singularity exec $SINGULARITY_IMAGE reconstructPar -latestTime | tee log.reconstructPar
singularity exec $SINGULARITY_IMAGE postChannel -latestTime | tee log.postChannel

#--- Final commands
echo "OpenFOAM script has reached the end"

```
{: .source}

> ## Important Part 1:
>
> ```bash
> #!/bin/bash --login
>
> #SBATCH --job-name=mpi-openfoam-training
> #SBATCH --partition=work
> #SBATCH --reservation=ContainersTraining
> #SBATCH --nodes=1
> #SBATCH --ntasks=8
> #SBATCH --ntasks-per-node=8
> #SBATCH --cpus-per-task=1
> #SBATCH --time=00:05:00
>
> #--- Load the singularity module (Pawsey's mpi-settings flavour):
> module load singularity/4.1.0-mpi
>
> #--- Using user's own image:
> export SINGULARITY_IMAGE="$MYSOFTWARE/singularity/images/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif"  #Adapt path and name to the correct ones
> echo "Using openfoam singularity image:"
> echo "SINGULARITY_IMAGE=$SINGULARITY_IMAGE"
> ```
> {: .source}
- The job reserves resources for 8 MPI tasks
- The `singularity/4.1.0-mpi` module is loaded
- The image to be used is defined in our user defined variable `SINGULARITY_IMAGE`
{: .solution}

> ## Important Part 2:
>
> ```bash
> #--- Execute pre-processing tools:
> #(These pre-processing tools are serial by design)
> singularity exec $SINGULARITY_IMAGE blockMesh | tee log.blockMesh
> singularity exec $SINGULARITY_IMAGE renumberMesh -overwrite -constant | tee log.renumberMesh
> singularity exec $SINGULARITY_IMAGE decomposePar -cellDist -force | tee log.decomposePar
>
> ...
>
> #--- Execute post-processing tools:
> #(These post-processing tools are serial by design)
> singularity exec $SINGULARITY_IMAGE reconstructPar -latestTime | tee log.reconstructPar
> singularity exec $SINGULARITY_IMAGE postChannel -latestTime | tee log.postChannel
> ```
> {: .source}
- All the required serial tools for pre- and post-processing are called with `singularity exec $SINGULARITY_IMAGE ...`
- Note that the `tee` commands are running in the host, so there pipe works fine passing the output from the singularity exectution and there's no need to use the `bash -c ` trick
{: .solution}

> ## Important Part 3:
>
> ```bash
> #--- Execute the parallel solver:
> #(Solvers use MPI parallelism by design)
> srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c 1 \
>   singularity exec $SINGULARITY_IMAGE pimpleFoam -parallel | tee log.pimpleFoam
> ```
> {: .source}
- `srun` launches the 8 MPI tasks in Setonix interconnect
- The parallel solver `pimpleFoam` is called using `singularity exec $SINGULARITY_IMAGE pimpleFoam -parallel`
{: .solution}

### How does Singularity interplay with the MPI launcher?

We'll comment on the environment variable definitions soon, now let's focus on the set of commands that make the simulation happen.

In particular, the fourth command is the only one using multiple processors through MPI:

```
mpirun -n $NTASKS \
  singularity exec openfoam_v2012.sif \
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
export SINGULARITYENV_LD_LIBRARY_PATH="$MPICH_ROOT/lib:\$LD_LIBRARY_PATH"
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
  singularity exec openfoam_v2012.sif \
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
