---
title: "Other notable container tools"
teaching: 10
exercises: 0
questions:
objectives:
- Get an overview of other tools of interest for containers on HPC
keypoints:
- HPCCM can be useful to write image recipes that are compatible both with Docker and Singularity
- Other than Singularity, at the moment Sarus seems to only other active project in terms of container engines for HPC
---


### HPC Container Maker

[HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker) (HPCCM) is a Python tool developed by Nvidia, whose only purpose is to write recipe files for containers.

The most interesting aspect is that recipes are written in a container engine agnostic format, and then HPCCM can translate them both in Dockerfiles and Singularity def files, allowing to simply produce container images in both formats.

Another interesting feature is that HPCCM ships a set of so called `building blocks`, *i.e.* ready to use units that install common packages in containers in an optimised way. For instance, these include compilers, MPI libraries, scientific libraries and a few applications.

To give an idea of how a HPCCM recipe looks like, here is one for the `lolcow` example we have explored when building Singularity and Docker containers:

```
Stage0 += baseimage(image='ubuntu:18.04')

Stage0 += packages(ospackages=['fortune', 'cowsay', 'lolcat'])
Stage0 += environment(variables={'PATH': '/usr/games:$PATH'})

Stage0 += label(metadata={'Author': '"Pawsey Supercomputing Centre"'})
Stage0 += label(metadata={'Version': 'v0.0.1'})
```
{: .source}

You can cd into the demo directory:

```
$ cd $SC19/demos/11_lolcow_hpccm
```
{: .bash}

And then create the def file:

```
$ hpccm --recipe lolcow.py --format singularity
```
{: .bash}

and the Dockerfile:

```
$ hpccm --recipe lolcow.py --format docker
```
{: .bash}

Note how this recipe does not produce exactly the def file and Dockerfile we have discussed in previous episodes. It is just meant to convey a general idea of the working principle of HPCCM.

More information on HPCCM can be found in the [HPCCM docs](https://github.com/NVIDIA/hpc-container-maker/tree/master/docs).


### Podman

[Podman](https://podman.io) is an open-source container engine maintained by Red Hat. It has quite similar features to Docker, with two important differences:

* runs daemon-less, making it simpler to deploy;
* can be used in root-less mode (with some [limitations](https://github.com/containers/libpod/blob/master/rootless.md)), making it a bit friendlier for HPC.

Like Docker, it still doesn't provide support for MPI, schedulers, GPU.

One reason to consider it is for building container images in rootless mode, right within a HPC system. On the other hand, for HPC runtime it lacks relevant features when compared to Singularity, Shifter and Sarus. Notably, at the time of writing it is still a bit buggy.

Interestingly, the API is mostly identical to that of Docker, so that in principle one could just use `podman` by means of

```
$ alias docker=podman
```
{: .bash}


### NERSC Shifter

[NERSC Shifter](https://docs.nersc.gov/programming/shifter/overview/) is a container engine developed by NERSC for HPC. 

It complies with HPC security requirements by design, and features native support for MPI and schedulers; interestingly, it cannot run GPU applications. It cannot be used to build images, just to run them.

At the time of writing, its development seems a bit stalled, for which reason we suggest considering CSCS Shifter/Sarus instead.


### CSCS Shifter (or Shifter-NG)

[CSCS Shifter](https://user.cscs.ch/tools/containers/shifter/) is a fork of NERSC Shifter by CSCS. 

Most notably, it adds GPU support to the runtime engine. 

It is being deprecated by CSCS, as they have evolved the project into Sarus, see below.


### Sarus

[Sarus](https://user.cscs.ch/tools/containers/sarus/) is the latest incarnation of a container runtime engine by CSCS. 

It is fully compliant with the Docker image format (whereas it cannot run Singularity images), natively supports schedulers, MPI, and GPU applications. Then in terms of runtime features it is mostly equivalent to Singularity (although at the moment it doesn't offer a feature comparable to OverlayFS). However, to build container images, it relies on the users being able to run Docker somewhere else. Also, uptake at the moment is quite limited compared to Singularity.

As at the moment Sarus seems the only practical alternative to running containers on HPC, let us proceed with a quick overview of the syntax. The key commands are:
* `sarus pull`: download container images (image name specification same as in Docker);
* `sarus run`: execute commands in containers;
* `sarus images`: list downloaded images;
* `sarus rmi`: remove downloaded image.

If you want to test it, you might just use the image `ubuntu:18.04` as a test bed, similar to what we did earlier on with Singularity and Docker.

Note how, at the time of writing, Sarus is still in a phase of beta testing, with the code base still subject to relevant changes.
