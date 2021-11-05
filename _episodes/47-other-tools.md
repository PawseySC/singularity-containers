---
title: "Other notable container tools"
teaching: 15
exercises: 0
questions:
objectives:
- Get an overview of other tools of interest for building containers
- Get an overview of other tools of interest for running containers on HPC
keypoints:
- HPCCM can be useful to write image recipes that are compatible both with Docker and Singularity
- Spack can be useful to automate generation of image recipes, including also micro-architecture optimisations
- In addition to Singularity, other interesting container engines for HPC exist
---


### HPC Container Maker

[HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker) (HPCCM) is a Python tool developed by Nvidia, whose only purpose is to write recipe files for containers.

The most interesting aspect is that recipes are written in a container engine agnostic format, and then HPCCM can translate them both in Dockerfiles and Singularity def files, allowing to simply produce container images in both formats.

Another interesting feature is that HPCCM ships a set of so called *building blocks*, *i.e.* ready to use units that install common packages in containers in an optimised way.  For instance, these include compilers, MPI libraries, scientific libraries and a few applications.

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
$ cd $TUTO/demos/lolcow_hpccm
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

Note how this recipe does not produce exactly the def file and Dockerfile we have discussed in previous episodes.  It is just meant to convey a general idea of the working principle of HPCCM.

A price to pay in using HPCCM is that translating a bash installation into a HPCCM recipe requires a lot of rewrite, as opposed to the smaller edits when porting to Dockerfiles and Singularity def files.

More information on HPCCM can be found in the [HPCCM docs](https://github.com/NVIDIA/hpc-container-maker/tree/master/docs).


### Spack

[Spack](https://spack.io) is a package manager for HPC, with the main purpose of automating the build from source of scientific applications and their dependencies.  It has many features and functionalities, with a relatively concise command interface to select compilers, package versions, dependencies, build options, compiler optimisations, and more.

Interestingly, Spack can also be used to automate the generation of Dockerfiles and Singularity def files.  In this regard, compared to other tools it easily allows to optimise the container build for a given CPU micro-architecture.  Thus, it can enable to reproducibly generate collections of container images for a given application, enforcing both portability and performance.

Another advantageous feature is its registry of over 5000 package recipes that are ready for use.

To give an example of how to generate a Dockefile with Spack, let's consider the bioinformatics package BLAST; suppose we need version 2.9.0.  Also suppose we want to optimise the application for running on Intel Haswell CPUs (or newer).  What we need is a so called Spack Environment file, in YAML format, which has to be called `spack.yaml`:

```
spack:
  specs:
  - blast-plus@2.9.0 target="haswell"

  container:
    format: "docker"

    images:
      os: "ubuntu:18.04"
      spack: "0.16"

    os_packages:
      build:
      - cpio
      final:
      - libgomp1

    labels:
      maintainer: "Pawsey Supercomputing Centre"
```
{: .source}

Without willing to provide an exaustive explanation of this file, note how the request for BLAST 2.9.0 optimised for Haswell is stated: `blast-plus@2.9.0 target="haswell"`.  Also note that we're requesting a Dockerfile by means of `format: "docker"`.  Now, you can cd into the demo directory:

```
$ cd $TUTO/demos/spack_blast
```
{: .bash}

And then create the Dockerfile just with:

```
$ spack containerize >Dockerfile
```
{: .bash}

More information on how to customise the creation of Dockefiles and def files with Spack can be found at the [Spack docs on container images](https://spack.readthedocs.io/en/latest/containers.html).


### Podman

[Podman](https://podman.io) is an open-source container engine maintained by Red Hat.  It has quite similar features to Docker, with some important differences:

* runs daemon-less, making it simpler to deploy in a shared environment such as HPC;
* can be used in root-less mode (with some [limitations](https://github.com/containers/libpod/blob/master/rootless.md)), making it a bit friendlier for HPC, including for image building right within a HPC system;
* provides basic support for MPI applications, see this [Podman blog post](https://podman.io/blogs/2019/09/26/podman-in-hpc.html);
* provides beta support for GPU applications, see this [Red Hat blog post](https://www.redhat.com/en/blog/how-use-gpus-containers-bare-metal-rhel-8).

These features are making Podman an increasingly appealing solution for running containers in HPC.

Interestingly, the API is mostly identical to that of Docker, so that in principle one could just use `podman` by means of

```
$ alias docker=podman
```
{: .bash}


### Charliecloud

[Charliecloud](https://hpc.github.io/charliecloud) is a promising container engine developed by LANL for HPC.

At the time of writing the project is in active development, with new features being added, some of them still in an experimental form.  

Charliecloud supports both MPI and GPUs.  
One of the most appealing features is the active effort in moving towards a completely unprivileged mode of operation with containers.  If Docker is available on the system, Charliecloud will use it as a privileged backend for pulling and building images.  
However, if no Docker is found (always the case in HPC), Charliecloud will fallback to its own unprivileged backends for pulling (experimental) and building.  Once the implementation of this fully unprivileged workflow reaches a stable status, Charliecloud will allow to run the entire container lifecycle on a HPC system.

Charliecloud is not a single binary application, instead it offers a collection of several executables for various purposes.  The same task (pulling, building, ..) can in general be performed in different ways.  
Here we're just showing one possible sequence of commands to pull and use the Ubuntu image `ubuntu:18.04`:

* Pull Ubuntu with `ch-image pull`
    ```
    $ ch-image pull ubuntu:18.04
    ```
    {: .bash}

    ```
    initializing storage directory: v2 /var/tmp/ubuntu/ch-image
    pulling image:    ubuntu:18.04
    requesting arch:  amd64
    manifest list: downloading: 100%
    manifest: downloading: 100%
    config: downloading: 100%
    layer 1/1: 2840553: downloading: 25.5/25.5 MiB (100%)
    flattening image
    layer 1/1: 2840553: listing
    validating tarball members
    resolving whiteouts
    layer 1/1: 2840553: extracting
    image arch: amd64
    done
    ```
    {: .output}

* Run a simple command with `ch-run`
    ```
    $ ch-run /var/tmp/ubuntu/ch-image/img/ubuntu:18.04 cat /etc/os-release
    ```
    {: .bash}


### Shifter

[Shifter](https://docs.nersc.gov/development/shifter/) is a container engine developed by NERSC for HPC.

It complies with HPC security requirements by design, and features native support for MPI and schedulers; at the moment, it cannot run GPU applications.  It cannot be used to build images, just to run them.  At the time of writing, it is quite popular within the HPC community in the USA.

Here are the key commands to use Shifter:
* `shifterimg pull`: download container images (image name specification same as in Docker);
* `shifter --image=<IMAGE>`: execute commands in containers;
* `shifterimg images`: list downloaded images;
* at the moment, Shifter does not have a shell command to remove container images (this has to be done by system administrators).

If you have access to a cluster with Shifter, here are some basic tests with the image `ubuntu:18.04`:

* Pull Ubuntu
    ```
    shifterimg pull ubuntu:18.04
    ```
    {: .bash}

    ```
    2020-11-02T22:38:46 Pulling Image: docker:ubuntu:18.04, status: READY
    ```
    {: .output}

* Run a simple command
    ```bash
    shifter --image=ubuntu:18.04 cat /etc/os-release
    ```


<!--
### CSCS Shifter (or Shifter-NG)

[CSCS Shifter](https://user.cscs.ch/tools/containers/shifter/) is a fork of NERSC Shifter by CSCS.

Most notably, it adds GPU support to the runtime engine.

It is being deprecated by CSCS, as they have evolved the project into Sarus, see below.
-->


### Sarus

[Sarus](https://sarus.readthedocs.io) is a container engine developed for HPC by CSCS in Switzerland.  It started as a fork of Shifter.

It is fully compliant with the Docker image format (whereas it cannot run Singularity images), natively supports schedulers, MPI, and GPU applications.  Then in terms of runtime features it is mostly equivalent to Singularity (although at the moment it doesn't offer a feature comparable to OverlayFS).  However, to build container images, it relies on the users being able to run Docker somewhere else.  Also, uptake at the moment is quite limited compared to Singularity.

Let us have a quick overview of the Sarus syntax.  The key commands are:
* `sarus pull`: download container images (image name specification same as in Docker);
* `sarus run`: execute commands in containers;
* `sarus images`: list downloaded images;
* `sarus rmi`: remove downloaded image.

If you want to test it, you might just use the image `ubuntu:18.04` as a test bed, similar to what we did earlier on with Singularity and Docker:

* Pull Ubuntu 
    ```
    $ sarus pull ubuntu/18.04
    ```
    {: .bash}

    ```
    # image            : index.docker.io/library/ubuntu:18.04
    # cache directory  : "/home/ubuntu/.sarus/cache"
    # temp directory   : "/tmp"
    # images directory : "/home/ubuntu/.sarus/images"
    > save image layers ...
    
    [..]
    
    > expanding image layers ...
    
    [..]
    
    > make squashfs image: "/home/ubuntu/.sarus/images/index.docker.io/library/ubuntu/18.04.squashfs"
    ```
    {: .output}

* Run a simple command
    ```
    $ sarus run ubuntu/18.04 cat /etc/os-release
    ```
    {: .bash}


### Enroot

[Enroot](https://github.com/NVIDIA/enroot) is Nvidia way of deploying containerised applications on their platforms.  

The tool has been kept simple by design, and can run fully unprivileged.  
It comes with native support for GPUs and Mellanox Infiniband (unsurprisingly).  Support for the Slurm scheduler is provided through the [Pyxis](https://github.com/NVIDIA/pyxis) plugin.  However, there is currently no mention of any MPI integration in the documentation.

It's definitely worth watching Nvidia's space for updates on this promising project.
