---
title: "Get ready on Setonix"
teaching: 10
exercises: 5
questions:
- "How can I connect to Setonix from Linux, macOS, or Windows?"
- "How can I request an interactive allocation for the training?"
- "How can I clone the training repository into my scratch directory?"
- "How can I start downloading container images needed later in the training?"
objectives:
- "Connect to a Setonix login node using SSH"
- "Request a compute node through an interactive Slurm allocation"
- "Clone the training repository into the Setonix scratch filesystem"
- "Set the environment variable used to locate the training files"
- "Launch Slurm jobs that pull container images needed later in the training"
keypoints:
- "Connect to Setonix using SSH and your Pawsey account"
- "Use login nodes for access and lightweight preparation, not computational work"
- "Run the hands-on exercises on a compute node obtained through Slurm"
- "Keep the training repository in your scratch directory"
- "Store reusable Singularity images in your personal image library"
---

#### Connect to Setonix

Before starting the hands-on exercises, connect to a Setonix login node. For this training, you should have been provided with credentials to make use of the resources reservation allocated for this course. If you are reading this outside of the training session, you would then need a Pawsey account that belongs to an active Setonix project allocation.

Setonix is accessed through Secure Shell (SSH) using the hostname `setonix.pawsey.org.au`. Replace `<username>` in the commands below with your Pawsey username and do not include the angle brackets.

The first time you connect, SSH may ask whether you trust the host key. Check that the hostname is `setonix.pawsey.org.au`, then enter `yes` to continue. When prompted, enter your Pawsey password. The password is not displayed while you type it.

##### Linux

Open a terminal and run:

```bash
$ ssh <username>@setonix.pawsey.org.au
```
{: .source}

##### macOS

Open _Terminal_, located in _Applications → Utilities_, and run:

```bash
$ ssh <username>@setonix.pawsey.org.au
```
{: .source}

##### Windows

Open _PowerShell_. In the Windows instructions, `PS>` represents the PowerShell prompt and must not be typed as part of the command.

Run:

```powershell
PS> ssh <username>@setonix.pawsey.org.au
```
{: .source}

Current versions of Windows normally include the OpenSSH client. If PowerShell reports that `ssh` is not recognised, enable the _OpenSSH Client_ optional feature in Windows or use an SSH client approved for your computer.

After a successful connection, you will be working on one of the Setonix login nodes. You can confirm the system hostname with:

```bash
$ hostname
```
{: .source}

```text
setonix-01
```
{: .output}

To end the SSH session when you have finished, run:

```bash
$ exit
```
{: .source}

> ## Login nodes are shared resources
>
> Use the login node to edit files, transfer data, clone repositories, and submit or monitor jobs. Do not run the computational exercises directly on a login node. Request a compute node through Slurm before running the hands-on workload.
{: .callout}

For additional connection details and troubleshooting, see Pawsey's [How to log into Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51926034/How+to+log+into+Setonix) documentation.

#### Request an interactive allocation

If you're running this tutorial on a shared system (_e.g._ Setonix at Pawsey), you should use one of the compute nodes rather than the login node. You can do this by requesting an interactive allocation from the scheduler, for instance on Setonix with Slurm (do this if you are not in an `salloc` interactive session yet):

```bash
$ salloc -N 1 -n 1 -c 16 --reservation=ContainersTraining -t 4:00:00
salloc: Granted job allocation 3453895
salloc: Waiting for resource configuration
salloc: Nodes nid000152 are ready for job
```
{: .source}

The job allocation number and compute-node hostname in your output will be different. Keep this terminal open while completing the hands-on exercises. Exiting the shell ends the interactive allocation.

#### Get ready for the hands-on

Before we start, let us ensure we have got the required files to run the tutorials.

If you haven't done so already, move to a suitable working directory and download the following GitHub repository. On Pawsey systems, use your scratch directory; on other HPC or cloud systems, use the equivalent working directory recommended by the system administrators.

```bash
$ cd "$MYSCRATCH"    # On Pawsey systems
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
```
{: .source}

Now move to the `demos` directory:

```bash
$ cd "${TUTO}/demos"
$ pwd
```
{: .source}

The working directory should be something like:

```output
/path/to/your/scratch/singularity-containers/demos
```
{: .output}

##### Start downloading images used later

Some container images used later in the training take time to download. Start their pulling jobs now so that the downloads can proceed through Slurm while the introduction is being taught.

Define the personal image library used throughout this training and create it if it does not already exist:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```
{: .source}

The launcher in the `demos` directory submits one independent Slurm job for each image. It waits two minutes between submissions so that, when training resources are limited, other participants have an opportunity to start their first download before one participant submits the next one.

Run the launcher:

```bash
$ ./launch_image_pulls.sh
```
{: .source}

The launcher initially requests the following images:

- the OpenFOAM image from Pawsey's Quay registry, saved as `${MY_LOCAL_LIBRARY}/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif`
- the Trinity image from Docker Hub, saved as `${MY_LOCAL_LIBRARY}/trinityrnaseq--2.8.6.sif`

You do not need to wait for the downloads to finish before continuing. Check the jobs with:

```bash
$ squeue --me
```
{: .source}

Each job writes its messages to a file named `pull-image-<job-id>.out`. Replace `<job-id>` with the Slurm job number reported when the job is submitted.
