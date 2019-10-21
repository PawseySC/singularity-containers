---
title: "Build your own container image with Docker"
teaching: 20
exercises: 10
questions:
objectives:
- Learn what is a Dockerfile and its basic syntax
- Learn how to build a container and push it to a web registry
keypoints:
- "A Dockerfile is a recipe that uses specific instructions to direct the image building process"
- "`docker build` is used to build images"
- "`docker push` is used to push images to a web registry"
---

### What is a Dockerfile? ###

A Dockerfile is a recipe to build a container image with Docker. It is basically a collection of the standard shell commands you would use to build your software through prompt; in addition, it contains Docker-specific instructions that handle the build process. We will see some examples below.


### Let's write a Dockerfile ###

We will build a very simple container image: a Ubuntu box featuring tools for building software and a text editor. Its Dockerfile will contain most of the basic Docker instructions that can also be used to build more complicated images.

First let us create a directory where we'll store the Dockerfile. This directory will be the so called Docker **build context**. Docker will include files in this directory in the build process and in the final image. As a by-product, this will make the build process longer and the image larger, so that we want to include only those strictly required for the build, even none if possible.

```
$ mkdir build_dockerfile
$ cd build_dockerfile
```
{: .bash}

Now use your favourite text editor to create a file named `Dockerfile` and edit it. Here is its contents:

```
FROM ubuntu:18.04
  
MAINTAINER Your Name <youremail@yourdomain>

RUN apt-get update && \
    apt-get install -y \
        autoconf \
        automake \
        g++ \
        gcc \
        gfortran \
        make \
        nano \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

VOLUME ["/data"]
WORKDIR /data

CMD ["/bin/bash"]
```
{: .source}

* `FROM`: compulsory, it provides the starting image we will use to build our customised one;
* `MAINTAINER`: details of the person who wrote the Dockerfile, optional;
* `RUN`: this is the most used instruction, that alllows to run most shell commands during the build. Multiple `RUN` instructions are often found in a single Dockerfile;
* `VOLUME`: creates a mount point ready to be used for mounting external (e.g. host) volumes; creates the corresponding directory if not existing;
* `WORKDIR`: changes directory to the specified path; the last current directory in the build will be the working directory in the running container.  
  **Note**: if you use instead `RUN cd <..>`, the changed directory will only persist within that `RUN` instruction, and then be lost in subsequent build steps;
* `CMD`: specifies the default command to be executed with the container. `bash` is the default anyway for Ubuntu containers, but it's good to be aware of this syntax.


### Building the image ###

Once the Dockerfile is ready, let us build the image with `docker build`:

```
$ docker build -t build:2019May29 .
```
{: .bash}

```
Sending build context to Docker daemon  2.048kB
Step 1/6 : FROM ubuntu:18.04
 ---> d131e0fa2585
[..]
Step 6/6 : CMD ["/bin/bash"]
 ---> Running in fb003b87b020
Removing intermediate container fb003b87b020
 ---> 8986ee76d9a9
Successfully built 8986ee76d9a9
Successfully tagged build:2019May29
```
{: .output}

In the command above, `.` is the location of the build context (i.e. the directory for the Dockerfile).  
The `-t` flag is used to specify the image name (compulsory) and tag (optional).

Any lowercase alphanumeric string can be used as image name; here we've used `build`. The image tag (following the colon) can be optionally used to maintain a set of different image versions on Docker Hub, and is a key feature in enabling reproducibility of your computations through containers; here we've used `2019May29`.

Adding the prefix `<Your Docker Hub account>/` to the image name is also optional and allows to push the built image to your Docker Hub registry (see below). 

The complete format for the image name looks like: `<Your Docker Hub account ^>/<Image name>:<Image tag ^>`. `^`These are optional.
 

### Layers in a container image ###

Note how the `RUN` instruction above is used to execute a sequence of commands to:

* update the list of available packages
* install a set of Linux packages
* clean build directories

We have concatenated all these commands in one using the `&&` linux operator, and then the `\` symbol to break them into multiple lines for readability.

We could have used one `RUN` instruction per command, so why concatenating instead?

Well, each `RUN` creates a distinct **layer** in the final image, increasing its size. It is a good practice to use as few layers, and thus `RUN` instructions, as possible, to keep the image size smaller.


### More Dockerfile instructions ###

Several other instructions are available, that we haven't covered in this introduction. You can find more information on them at [Dockerfile reference](https://docs.docker.com/engine/reference/builder/). Just to mention a few possibilities:

* `ARG`: set temporary values that will be used during the build process, and that might need to be changed in future builds; a common use is to specify package versions; `docker build` has an option to change at build time the value of temporary `ARG` variables set in the Dockerfile: `--build-arg <variable>=<value>`;
* `ENV`: set environment variables that will persist at runtime in the container; **DO NOT** use `RUN export <..>` to this end, as the variable will be lost after the `RUN` step is completed;
* `ADD`/`COPY`: embed files/directories from your computer into the container image;
* `EXPOSE`: make the container listen on specified network ports;
* `CMD`/`ENTRYPOINT`: tweak the default behaviour of the executing container;
* `USER`: switch user.


### Pushing the image to Docker Hub ###

If you have a (free) Docker Hub account you must first login to Docker.

```
$ docker login
```
{: .bash}


You are now ready to push your newly created image to the Docker Hub web registry.

First, let us create a second tag for the image, that includes your Docker Account. To this end we'll use `docker tag`:

```
$ docker tag build:2019May29 <your-dockerhub-account>/build:2019May29
```
{: .bash}

Now we can push the image:

```
$ docker push <your-dockerhub-account>/build:2019May29
```
{: .bash}

```
The push refers to repository [docker.io/marcodelapierre/build]
cab15c00fd34: Pushed 
cf5522ba3624: Pushed 
[..]
2019May29: digest: sha256:bcb0e09927291c7a36a37ee686aa870939ab6c2cee2ef06ae4e742dba4bb1dd4 size: 1569
```
{: .output}

Congratulations! Your image is now publicly available for anyone to pull.


### Base images for Python ###

[continuumio/miniconda2](https://hub.docker.com/r/continuumio/miniconda2/tags) and [continuumio/miniconda3](https://hub.docker.com/r/continuumio/miniconda3/tags) are Docker images provided by the maintainers of the [Anaconda](https://anaconda.org) project. They ship with Python 2 and 3, respectively, as well as `pip` and `conda` to install and manage packages. At the time of writing, the most recent version is `4.5.12`, which is based on Python `2.7.15` and `3.7.1`, respectively.

Among other use cases, these base images can be very useful for maintaining Python containers, as well as bioinformatics containers based on the [Bioconda](https://bioconda.github.io) project.

If you need interactive Jupyter Notebooks, [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/) are a series of dedicated container images. Among others, there is the base SciPy image [jupyter/scipy-notebook](https://hub.docker.com/r/jupyter/scipy-notebook/tags/), the data science image [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook/tags/), and the machine learning image [jupyter/tensorflow-notebook](https://hub.docker.com/r/jupyter/tensorflow-notebook/tags/).


### Base images for R ###

The [Rocker Project](https://www.rocker-project.org) maintains a number of good R base images. Of particular relevance is [rocker/tidyverse](https://hub.docker.com/r/rocker/tidyverse/tags), which embeds the basic R distribution, an RStudio web-server installation and the [tydiverse](https://www.tidyverse.org) collection of packages for data science, that are also quite popular across the bioinformatics community of [Bioconductor](https://www.bioconductor.org). At the time of writing, the most recent version is `3.5.3`.

Other more basic images are [rocker/r-ver](https://hub.docker.com/r/rocker/r-ver/tags) (R only) and [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio/tags) (R + RStudio).


> ## Build your own Scientific Python container ##
> 
> Using `continuumio/miniconda3:4.5.12` as base image, create an image called `mypython`, which includes the Python packages **numpy**, **scipy** and **pandas**. Hint: you can use `pip install` or `conda install -y` in the Dockerfile to this end. 
> 
> If you have a Docker Hub account, for the image name use the format `<Your Docker Hub account>/<Image name>:<Version tag>`. Then, push the image to the web registry.
> 
> > ## Solution ##
> > 
> > Dockerfile:
> > 
> > ```
> > FROM continuumio/miniconda3:4.5.12
> > 
> > MAINTAINER Marco De La Pierre <marco.delapierre@pawsey.org.au>
> > 
> > RUN pip install numpy scipy pandas
> > 
> > VOLUME ["/data"]
> > WORKDIR /data
> > 
> > CMD ["/bin/bash"]
> > ```
> > {: .source}
> > 
> > Build the image:
> > 
> > a) Plain (no Docker Hub account):
> > 
> > ```
> > $ docker build -t mypython:2019Apr23 .
> > ```
> > {: .bash}
> > 
> > b) With a Docker Hub account:
> > 
> > ```
> > $ docker build -t <your-dockerhub-account>/mypython:2019Apr23 .
> > ```
> > {: .bash}
> > 
> > Push the image (optional):
> > 
> > ```
> > $ docker push <your-dockerhub-account>/mypython:2019Apr23
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


> ## Best practices ##
> 
> * for stand-alone packages, it is suggested to use the policy of one container per package
> * for Python or R pipelines, it may be handier to use the policy of a single container for the entire pipeline
> * [Best practices for writing Dockerfiles](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/) are found in the Docker website
{: .callout}
