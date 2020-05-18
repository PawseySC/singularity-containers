---
title: "Building images with Docker"
teaching: 15
exercises: 5
questions:
objectives:
- Learn the pros&cons of building with Singularity vs Docker
- Learn how to build and share an image with Docker
keypoints:
- Build images starting from a `Dockerfile` recipe with `docker build`
- Push images to a web registry with `docker push`
---


### Why using Docker for container builds?

Over the whole of this tutorial we're proposing Singularity as the principal tool to handle and run containerised applications in HPC.  Do we really need to extend our toolkit to include [Docker](https://hub.docker.com/search/?type=edition&offering=community)?

To better inform an answer to this question, here are some of the advantages when building with one or the other tool.

#### Singularity
* Single file image, can be handled as any other file and shared easily
* Unambiguous container usage modes, via distinct keywords: `exec`, `shell`, `run`, `instance` (see episode on GUI applications)
* Ability to sign/verify images for improved security

#### Docker
* Compatibility: image format can be run by all existing container engines
* Quick development: layered image format allows caching, for reduced build time during repeated builds

Note how, at present, both tools require root privileges for building, implying that this step cannot be performed on HPC, and requires a dedicated machine instead.

Although Singularity builds offer some interesting advantages, there's a single item that right now makes using Docker preferred in most situations.  
It's **compatibility**.  Build with Docker, and you'll know the resulting image can be run from every container engine anywhere in the world.


### A Dockerfile recipe

Let's cd into the relevant demo directory:

```
$ cd $TUTO/demos/11_lolcow_docker
```
{: .bash}

It can be interesting to have an idea of how to build images with Docker. In fact, as we mentioned earlier on, the layered image format of Docker can sometimes help in reducing image development time. In addition, Docker images are quite universally compatible, as they can be run by Singularity, too.

We're going to build a very similar image to the one we built with Singularity. The `Dockerfile` recipe file in the demo directory looks like:

```
FROM ubuntu:18.04

MAINTAINER Pawsey Supercomputing Centre

RUN apt-get -y update && \
  apt-get -y install fortune cowsay lolcat

ENV PATH=/usr/games:$PATH

VOLUME /data
WORKDIR /data

CMD fortune | cowsay | lolcat
```
{: .source}

The directory where the the Dockerfile is stored is the so called the Docker **build context**. Docker will include files in this directory in the build process and in the final image. As a by-product, this will make the build process longer and the image larger, so that we want to include only those strictly required for the build, even none if possible.

Let's comment on the Docker instructions that appear in this Dockerfile.

* `FROM`: compulsory, it provides the starting image we will use to build our customised one;
* `MAINTAINER`: details of the person who wrote the Dockerfile, optional;
* `RUN`: this is the most used instruction, that allows to run most shell commands during the build. Multiple `RUN` instructions are often found in a single Dockerfile;
* `ENV`: set environment variables that will persist at runtime in the container; **DO NOT** use `RUN export <..>` to this end, as the variable will be lost after the `RUN` step is completed;
* `VOLUME`: creates a mount point ready to be used for mounting external (e.g. host) volumes; creates the corresponding directory if not existing;
* `WORKDIR`: changes directory to the specified path; the last current directory in the build will be the working directory in the running container.  
  **Note**: if you use instead `RUN cd <..>`, the changed directory will only persist within that `RUN` instruction, and then be lost in subsequent build steps;
* `CMD`: specifies the default command to be executed with the container, in case no other command is provided.

More information on the Dockerfile syntax can be found at the [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).


### Layers in a container image

Note how the `RUN` instruction above is used to execute a sequence of commands to update the list of available packages and install a set of Linux packages.  
We have concatenated all these commands in one using the `&&` linux operator, and then the `\` symbol to break them into multiple lines for readability.

We could have used one `RUN` instruction per command, so why concatenating instead?  
Well, each `RUN` creates a distinct **layer** in the final image, increasing its size. It is a good practice to use as few layers, and thus `RUN` instructions, as possible, to keep the image size smaller.


### Building the container image

Once the Dockerfile is ready, let us build the image with `docker build`:

```
$ sudo docker build -t lolcow:1Nov19 .
```
{: .bash}

In the command above, `.` is the location of the build context (i.e. the directory for the Dockerfile).  
The `-t` flag is used to specify the image name (compulsory) and tag (optional).

Any lowercase alphanumeric string can be used as image name; here we've used `lolcow`. The image tag (following the colon) can be optionally used to maintain a set of different image versions on Docker Hub, and is a key feature in enabling reproducibility of your computations through containers; here we've used `1Nov19`.

Adding the prefix `<Your Docker Hub account>/` to the image name is also optional and allows to push the built image to your Docker Hub registry (see below). 

The complete format for the image name looks like: `<Your Docker Hub account ^>/<Image name>:<Image tag ^>`. `^`These are optional.

This is the output of our build:

```
Sending build context to Docker daemon  2.048kB
Step 1/7 : FROM ubuntu:18.04
 ---> 775349758637
Step 2/7 : MAINTAINER Pawsey Supercomputing Centre
 ---> Running in 91c109dfd5ba
Removing intermediate container 91c109dfd5ba
 ---> 361490204a2c
Step 3/7 : RUN apt-get -y update &&   apt-get -y install fortune cowsay lolcat
 ---> Running in 4543b6bb99f1

[..]

Removing intermediate container 4543b6bb99f1
 ---> 7958a569068f
Step 4/7 : ENV PATH=/usr/games:$PATH
 ---> Running in 86282799c41f
Removing intermediate container 86282799c41f
 ---> 3ffdfe179e34
Step 5/7 : VOLUME /data
 ---> Running in f93de5446caa
Removing intermediate container f93de5446caa
 ---> 9c174e36bf3a
Step 6/7 : WORKDIR /data
 ---> Running in eed67d591239
Removing intermediate container eed67d591239
 ---> 36cc09b2c59b
Step 7/7 : CMD fortune | cowsay | lolcat
 ---> Running in 87a464d2ee67
Removing intermediate container 87a464d2ee67
 ---> 3c62a0f2e06e
Successfully built 3c62a0f2e06e
Successfully tagged lolcow:1Nov19
```
{: .output}

Let's give this image a go! Let's execute it without any argument to use the default command:

```
$ sudo docker run --rm lolcow:1Nov19
```
{: .bash}

```
 _______________________________________
/ Good news. Ten weeks from Friday will \
\ be a pretty good day.                 /
 ---------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}

Note how the default command can be readily overwritten:

```
$ sudo docker run --rm lolcow:1Nov19 echo "Hello World!"
```
{: .bash}

```
Hello World!
```
{: .output}


### Pushing the image to Docker Hub

If you have a (free) Docker Hub account you must first login to Docker.

```
$ sudo docker login
```
{: .bash}


You are now ready to push your newly created image to the Docker Hub web registry.

First, let us create a second tag for the image, that includes your Docker Account. To this end we'll use `docker tag`:

```
$ sudo docker tag lolcow:1Nov19 <your-dockerhub-account>/lolcow:1Nov19
```
{: .bash}

Now we can push the image:

```
$ sudo docker push <your-dockerhub-account>/lolcow:1Nov19
```
{: .bash}

```
The push refers to repository [docker.io/marcodelapierre/lolcow]
9d2959e72647: Pushed 
317d47a452af: Pushed 
e0b3afb09dc3: Mounted from library/ubuntu 
6c01b5a53aac: Mounted from library/ubuntu 
2c6ac8e5063e: Mounted from library/ubuntu 
cc967c529ced: Mounted from library/ubuntu 
1Nov19: digest: sha256:295c5695e2b05f6123bc2d8669ec7b66e45df5000ab9fc45ce3566ae3c0d839e size: 1571
```
{: .output}

Your image is now publicly available for anyone to pull.
