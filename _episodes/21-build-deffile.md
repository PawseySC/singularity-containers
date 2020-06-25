---
title: "Inside a build definition file"
teaching: 15
exercises: 5
questions:
objectives:
- Discuss what a definition file is
- Learn the basic syntax of a definition file
- Discuss some useful base container images
keypoints:
- Use `Bootstrap` and `From` to specify the build starting point
- The `%post` section contains the list of commands needed to install and setup packages in the image
- The `%environment` section allows you to specify shell variable definitions that are required at runtime
---


### What is a def file?

A *definition file*, or *def file*, is a recipe to build a container image with Singularity.  
In practice, it is a collection of the standard shell commands you would use to build your software, plus Singularity-specific header lines that handle the build process.  We will discuss these below with an example.  
Although there is no mandatory naming convention for def files, they are often characterised by the suffix `.def`.


### Core syntax of a definition file

We're going to have a closer look at the def file that was used in the introductory episode on building images.  If you've already gone through that episode, you'll be able to use the newly created image to test some new Singularity commands.  
This example is adapted from this well crafted [Singularity Tutorial](https://github.com/ArangoGutierrez/Singularity-tutorial).

Let us cd into the appropriate directory:

```
$ cd $TUTO/demos/lolcow
```
{: .bash}

How does `lolcow.def` look like?

```
Bootstrap: docker
From: ubuntu:18.04

%post
    apt-get -y update
    apt-get -y install fortune cowsay lolcat

%environment
    export PATH=/usr/games:$PATH

%labels
    Author Pawsey Supercomputing Centre
    Version v0.0.1

%help
    This is a demo container used to illustrate a def file.

%runscript
    fortune | cowsay | lolcat
```
{: .source}


The first line in the *def file* is `Bootstrap: docker`.  
This tells Singularity how the image has to be initialised. `docker` means that we are going to start with a base image from Docker Hub.  Another common way to bootstrap is using `library`, which will grab an image from the Sylabs Cloud.  The image is specified in the next line, in this case `From: ubuntu:18.04`.  
Note how we started from Ubuntu 18.04 in Docker Hub, not Sylabs Cloud, as the former version has got a bit of a richer, more useful configuration.

Next is a section that starts with the header `%post`.  This is basically a sequence of commands to be executed to install packages in the image, in essence the same commands you would use for installation in a Linux box.  Here we are ensuring we have an up-to-date list of packages, and then we are installing three Linux utilities.

The section `%environment` sets up environment variables that need to be defined at runtime rather than at build time.  Here the `PATH` needs to be updated to reflect the location of the three utilities that we installed in the `%post` section.

Although we are not using it in this *def file*, another section that is often useful can be defined by the header `%files`, like in:

```
%files
    <src-file> <dst-file>
```
{: .source}

This section is used to copy files from the host, *i.e.* <src-file>, inside the container in the destination <dst-file>.


### Documenting the container image

The `%labels` section is used to add metadata to the container image.  These can be then inspected by using

```
$ singularity inspect lolcow.sif
```
{: .bash}

```
==labels==
org.label-schema.usage.singularity.deffile.from: ubuntu:18.04
org.label-schema.usage.singularity.version: 3.3.0
Version: v0.0.1
org.label-schema.build-date: Tuesday_29_October_2019_14:44:19_UTC
org.label-schema.usage: /.singularity.d/runscript.help
org.label-schema.usage.singularity.runscript.help: /.singularity.d/runscript.help
Author: Pawsey Supercomputing Centre
org.label-schema.schema-version: 1.0
org.label-schema.usage.singularity.deffile.bootstrap: docker
```
{: .output}

See how the `Author` and `Version` metadata are in the list.

The text content of the `%help` section is also embedded in the image, and can be accessed via

```
$ singularity run-help lolcow.sif
```
{: .bash}

```
    This is a demo container used to illustrate a def file.


```
{: .output}

This can be useful to provide a description of the container, or even instructions on how to use it.

Finally, note how the def file used to generate the image can be displayed using

```
$ singularity inspect -d lolcow.sif
```
{: .bash}


### Run a container as an application

There's one section of the def file we haven't commented on yet.  `%runscript` allows you to define a default command for the image. This command can then be used if you run the container as an executable:

```
$ ./lolcow.sif
```
{: .bash}

```
 ______________________________________
/ A few hours grace before the madness \
\ begins again.                        /
 --------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}

Or, if you need to specify Singularity runtime flags, *e.g.*:

```
$ singularity run -B $TUTO/_episodes lolcow.sif
```
{: .bash}

```
 ___________________________________
< You will outgrow your usefulness. >
 -----------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}


> ## Use `%runscript` only for the talking cow!
> 
> Beside this amusing example, the use of `%runscript` is well documented and advertised in the Singularity docs.  However, we are here suggesting to avoid using such functionality for production container images.
> 
> Why?  
> Well, Singularity has two main ways to execute commands:
> * `> singularity exec <image> <cmd>`
> * `singularity run <image>` (default command only)
> 
> Now, most software packages provide multiple executables, depending on specific functions to be used, or tasks to be run.  There are some ill defined aspects here.  
> First, what type of general tip can we suggest for picking the *default* command to be associated with `run`?  It's hard to tell.  
> Second, this would result in an inconsistent user experience, where one command in an image is executed using `run`, and all the others via `exec`.  
> 
> For these reasons, our recommended best practice is to always use the `exec` syntax in Singularity to execute commands from a container.
{: .callout}


### Advanced build options

The def file specification has a number of other interesting features.  To know more about them you can visit the [Sylabs docs on def files](https://sylabs.io/guides/3.5/user-guide/definition_files.html).

1. In the episode on GUI enabled web applications we'll see how to use `%startscript` to configure the behaviour of containers running in background as **long running web servers**.

2. If you are in a development phase, where you don't know yet what you will include in your final container image, you can start with a **sandbox** image.  This is a special type of image designed for development purposes, consisting not of a single file, but instead of a directory.  To create one, run something like:

    ```
    $ sudo singularity build --sandbox playbox/ docker://ubuntu:18.04
    ```
    {: .bash}

    Then to open it and play, run:

    ```
    $ sudo singularity shell --writable playbox/
    ```
    {: .bash}

    Do NOT use sandboxes for production, as their design is not reproducible by nature.  
    More information on sandbox images can be found at the [Sylabs docs on building images](https://sylabs.io/guides/3.5/user-guide/build_a_container.html#creating-writable-sandbox-directories).

3. One last notable feature is the ability to use PGP keys to **sign and verify** container images.  In this way, users of 3rd party containers can double check that the image they're running is bit-by-bit equivalent to the one that the author originally built, largely reducing the possibility of running containers infected by malware.  
    You can find more on this topic at the [Sylabs docs on signing and verifying containers](https://sylabs.io/guides/3.5/user-guide/signNverify.html).


### Useful base images

At the time of writing, [Docker Hub](https://hub.docker.com) is the most popular web registry for general purpose container images.  Therefore all images mentioned below are hosted in this registry.

#### CUDA
[nvidia/cuda](https://hub.docker.com/r/nvidia/cuda) has images to build GPU enabled applications.  There are different image types for different needs.  Tags containing `runtime` are suitable for binary applications that are ready to run; if you need to compile GPU code, pick tags containing `devel` instead.  Different OS flavours are available, too.

#### MPI
As you can see in the episode on MPI applications, when containerising this type of software the MPI libraries in the image need to be ABI compatible with the MPI libraries in the host.  The Pawsey Supercomputing Centre maintains some dedicated base images at [pawsey/mpich-base](https://hub.docker.com/r/pawsey/mpich-base), for building images that will run on our HPC systems.

#### Python
[python](https://hub.docker.com/_/python) hosts the official Python images.  Different versions are available for some OS flavours.  Smaller base images have tags ending with `-slim`.

[continuumio/miniconda3](https://hub.docker.com/r/continuumio/miniconda3) are images provided by the maintainers of the [Anaconda](https://anaconda.org) project.  They ship with Python 3, as well as `pip` and `conda` to install and manage packages.

If you need interactive Jupyter Notebooks, [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/) maintain a series of dedicated container images.  Among others, there is the base SciPy image [jupyter/scipy-notebook](https://hub.docker.com/r/jupyter/scipy-notebook), the data science image [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook), and the machine learning image [jupyter/tensorflow-notebook](https://hub.docker.com/r/jupyter/tensorflow-notebook).

#### R
The [Rocker Project](https://www.rocker-project.org) maintains a number of good R base images.  Of particular relevance is [rocker/tidyverse](https://hub.docker.com/r/rocker/tidyverse), which embeds the basic R distribution, an RStudio web-server installation and the [tidyverse](https://www.tidyverse.org) collection of packages for data science.

Other more basic images are [rocker/r-ver](https://hub.docker.com/r/rocker/r-ver) (R only) and [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio) (R + RStudio).
