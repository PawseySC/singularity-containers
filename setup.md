---
layout: page
title: Setup
root: .
---

The main requirement for this workshop is a personal computer with a web browser and a command line shell program.  
These will allow you to follow the online materials and to login to a facility with the required software stack.


> ## SC19 Workshop Attendee Instructions
>
> You will will need username and password to access the training VMs, which will be provided by the instructors.
>
> You will find a green and pink post-it note at your spot.  On the green note is a 3 digit number ("XXX" below).  Your login details today follow a pattern of:
>
> **username**: userXXX  
> **password**: tutorialXXX
>
> You don't need to read further for today.
{: .callout}


### Regular users of this tutorial: read here

* **Zeus @Pawsey**: if you have access, Singularity can be loaded with `module load singularity`. MPI libraries are configured properly, and GPU applications can be run on the Slurm partition `gpuq`.

* **Nimbus Cloud @Pawsey**: if you have access, both Singularity and Docker are preinstalled in the *Ubuntu Pawsey* base image.
<!-- Test: Ubuntu 18.04 VM with 2 cores, 6 GB RAM, 40 GB disk -->

* **BYO Device**: if you have a Linux box, you can install the required software yourself (might take a while):

  * Essential (core of the tutorial)
    - Singularity : [script]({{ page.root }}/files/install-singularity.sh) \| [docs](https://sylabs.io/guides/3.3/user-guide/installation.html)

  * Desirable (to run all Singularity examples)
    - MPICH library : [script]({{ page.root }}/files/install-mpich.sh) \| [docs](https://www.mpich.org/documentation/guides/)
    - Nvidia GPU driver (GPU card required)
    - Slurm scheduler
    - Nextflow engine : [script]({{ page.root }}/files/install-nextflow.sh) \| [docs](https://www.nextflow.io/docs/latest/getstarted.html)

  * Optional (extra applications for last two episodes)
    - Docker : [docs (unofficial)](https://www.itzgeek.com/how-tos/linux/ubuntu-how-tos/how-to-install-docker-on-ubuntu-18-04-lts-bionic-beaver.html)
    - HPCCM : [script]({{ page.root }}/files/install-hpccm.sh) \| [docs](https://github.com/NVIDIA/hpc-container-maker/blob/master/docs/getting_started.md)
    - Podman : [script]({{ page.root }}/files/install-podman.sh) \| [docs](https://podman.io/getting-started/installation)
    - Sarus : [script]({{ page.root }}/files/install-sarus.sh) \| [docs](https://sarus.readthedocs.io/en/latest/install/requirements.html)

**Notes**
* Install scripts have been tested on a Ubuntu machine through a user that can run *sudo* commands without password prompts. There's no warranty they will work in your Linux box, you should consider them as templates.
* To install Singularity on a Mac or Windows machine, follow these instructions by Sylabs on [Setting up Singularity with Vagrant](https://sylabs.io/guides/3.3/user-guide/installation.html#install-on-windows-or-mac).
