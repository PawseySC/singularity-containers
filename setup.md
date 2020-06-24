---
layout: page
title: Setup
root: .
---


### Key requirement

The main requirement for this workshop is a personal computer with a web browser and a command line shell program.  
*Windows* users: get [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html), [Visual Studio Code](https://code.visualstudio.com/) or PuTTY.  
This equipment will allow you to follow the online materials and to login to a facility with the required software stack.

*Optional*: if you want to run the episode on X11 applications and you have a Windows box, you will also need to install [Cygwin/X](https://x.cygwin.com).  
If you have macOS, you will need [XQuartz](https://www.xquartz.org) installed instead.


### Self-paced: Pawsey users

* **Nimbus Cloud @Pawsey**: if you have access, both Singularity and Docker are preinstalled in the *Ubuntu Pawsey* base image.
<!-- Test: Ubuntu 18.04 VM with 2 cores, 6 GB RAM, 40 GB disk -->

* **Zeus @Pawsey**: if you have access, Singularity can be loaded with `module load singularity`. MPI libraries are configured properly, and GPU applications can be run on the Slurm partition `gpuq`.


> ## Advanced self-paced: installation DIY (Do It Yourself)
> 
> If you have your own machine, you can install the required software yourself (might take a while).  
> Note that you will need ***admin* privileges** in the machine to finalise the installation.
> 
> > ## Linux box: read here
> > 
> > * Essential (core of the tutorial)
> >   - Singularity : [script]({{ page.root }}/files/install-singularity.sh) \| [docs](https://sylabs.io/guides/3.5/user-guide/quick_start.html)
> > * Desirable (to run all the base episodes)
> >   - Docker : [script]({{ page.root }}/files/install-docker.sh) \| [docs (unofficial)](https://www.itzgeek.com/how-tos/linux/ubuntu-how-tos/how-to-install-docker-on-ubuntu-18-04-lts-bionic-beaver.html)
> >   - MPICH library : [script]({{ page.root }}/files/install-mpich.sh) \| [docs](https://www.mpich.org/documentation/guides/)
> >   - Nvidia GPU driver (GPU card required)
> >   - Slurm scheduler (can still run the MPI examples without it)
> > * Optional (bonus episodes)
> >   - Nextflow engine : [script]({{ page.root }}/files/install-nextflow.sh) \| [docs](https://www.nextflow.io/docs/latest/getstarted.html)
> >   - Environment modules : [script]({{ page.root }}/files/install-modules.sh) \| [docs](http://modules.sourceforge.net)
> >   - Docker Compose: [script]({{ page.root }}/files/install-dockercompose.sh) \| [docs](https://docs.docker.com/compose/)
> >   - HPCCM : [script]({{ page.root }}/files/install-hpccm.sh) \| [docs](https://github.com/NVIDIA/hpc-container-maker/blob/master/docs/getting_started.md)
> >   - Podman : [script]({{ page.root }}/files/install-podman.sh) \| [docs](https://podman.io/getting-started/installation)
> >   - Sarus : [script]({{ page.root }}/files/install-sarus.sh) \| [docs](https://sarus.readthedocs.io/en/latest/install/requirements.html)
> >   - Charliecloud : [script]({{ page.root }}/files/install-charliecloud.sh) (uses [Spack](https://spack.io)) \| [docs](https://hpc.github.io/charliecloud)
> >   - Enroot : [script]({{ page.root }}/files/install-enroot.sh) \| [docs](https://github.com/NVIDIA/enroot/blob/master/doc/installation.md)
> > 
> > **Note:** install scripts have been tested on a Ubuntu machine through a user that can run *sudo* commands without password prompts. There's no warranty they will work in your Linux box, you should consider them as templates.
> {: .solution}
> 
> #### macOS or Windows machine
> 
> For *Singularity*, you will need to setup a Linux virtual machine, and then follow the same instructions as above.  
> It's not as bad as it sounds... the main two options are:
>   - Vagrant: follow these instructions by Sylabs on [Setting up Singularity with Vagrant](https://sylabs.io/guides/3.5/admin-guide/installation.html#installation-on-windows-or-mac) (*macOS* users: DO NOT use the proposed *Singularity Desktop*, use Vagrant instead);
>   - Multipass: follow instructions from the [Multipass Homepage](https://multipass.run).
> 
> For *Docker*, you can download and run installers for [macOS](https://hub.docker.com/editions/community/docker-ce-desktop-mac/) and [Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows/).
> 
> #### All platforms
> 
> A more detailed discussion on how to setup Singularity can be found in a dedicated episode of this tutorial.
{: .challenge}
