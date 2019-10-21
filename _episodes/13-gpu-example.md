---
title: "Molecular dynamics with GPU containers"
teaching: 0
exercises: 20
questions:
objectives:
- Get started with Nvidia GPU containers for HPC applications
keypoints:
- You can use containers to ship GPU applications
---

### Nvidia GPU Cloud ###

**Note**: if your Docker machine has got an Nvidia GPU installed, then you can install the `nvidia-docker` (e.g. see <https://devblogs.nvidia.com/gpu-containers-runtime> ). For this tutorial, we are instead going to use Shifter and the GPU nodes available on the Zeus HPC system at Pawsey.

The GPU manufacturer Nvidia has a dedicated web registry for container images, shipping GPU optimised applications: <https://ngc.nvidia.com>.

To access this registry, you'll need a free account. Go to <https://ngc.nvidia.com>, complete the procedure to **Create an account**, then **Sign in** (currently both options are available on the top right corner of the page).

You can browse the available containerised packages by browsing the different categories. E.g. click on the **High Performance Computing** box, then click on the **Gromacs** one. The page will briefly discuss the code, with instructions on how to pull and run the container.

One key aspect is that to pull containers from the Nvidia GPU Cloud you'll need to provide Docker/Shifter with some login credentials. The username is always `$oauthtoken`, whereas the password needs to be retrieved in your Nvidia account. This password can be regenerated, so don't worry if you lose it. On the left side of the page, click on the **Configuration** tab, then on the button **Get API Key**. Then in the new page click on the button **Generate API Key** (currently at the top right corner), and click on confirm. An API Key string will appear in the page, copy it in your clipboard, then store it somewhere useful for your shell session, for instance in an environment variable:

```
$ export NVIDIA_KEY=<Paste Key Here>
```
{: .bash}

### Run a molecular dynamics simulation on GPUs with containers ### 

For our example we are going to use Gromacs, a quite popular molecuar dynamics package, among the ones that have been optimised to run on GPUs through Nvidia containers.

Let us pull the Gromacs container on the Zeus HPC system at Pawsey. We'll need to use `shifter pull --login` in order to enter the Nvidia cloud credentials. Also note we are prepending `nvcr.io/` to the repository name, to tell Shifter we are pulling from the Nvidia GPU Cloud:

```
$ module load shifter
$ sg $PAWSEY_PROJECT -c 'shifter pull --login nvcr.io/hpc/gromacs:2018.2'
```
{: .bash}

When prompted, enter `$oauthtoken` as username, and the copied key as password; if you've forgotten the latter, retrieve it from the environment variable you exported, using `echo $NVIDIA_KEY`.

Now, let us create a working directory, and pull some example files:

```
$ cd $MYSCRATCH
$ mkdir gpu_example
$ cd gpu_example

$ wget ftp://ftp.gromacs.org/pub/benchmarks/water_GMX50_bare.tar.gz
$ tar xzf water_GMX50_bare.tar.gz
$ cp water-cut1.0_GMX50_bare/1536/* .
```
{: .bash}

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

