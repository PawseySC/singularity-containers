---
title: "Overview of Docker"
teaching: 15
exercises: 15
questions:
objectives:
- Be aware of the pros and cons compared to Singularity
- Get started with Docker, download and run container images
- Discuss how to build and share an image with Docker
keypoints:
- Use `docker pull` to download container images
- Use `docker run` to execute commands in containers
- Build images starting from a `Dockerfile` recipe with `docker build`
- Push images to a web registry with `docker push`
---


### What's the deal with Docker?

[Docker](https://hub.docker.com/search/?type=edition&offering=community) has been the first container engine to get widespread popularity.  It has achieved this mostly in the world of IT companies, where it can be a very effective tool in the hands of system administrators, to deploy all sorts of micro-services.  It can also be a useful engine for running containers in laptops, personal workstations, and cloud VMs.  Among its advantages: 

* *root* execution allows for complete control and customisation;
* *isolation* over *integration*: by default Docker runs containers in complete isolation compared to the host, with highest security.  Users are in control of plugging in additional host features, such as directories/volumes, networks, communications ports;
* *docker-compose* to define and run multi-container applications, allowing to manage complex workflows; *e.g.* this can make Docker convenient for deploying long running services including Jupyter and RStudio servers;
* caching of exited containers, to eventually restart them;
* layered image format allows for caching of container building steps during build time, reducing development time.

On the other hand, some features make it not ideal for HPC.  These include:

* users need *root* privileges to run it, which is not really a good idea in a shared system;
* *isolation* over *integration* means users need to get used to a more articulated syntax to get things working with typical HPC applications;
* no support offered to interface Docker with MPI runtime, or HPC schedulers;
* usually requires an up-to-date kernel.

As you might encounter Docker in your container journey, let's have a quick look at how the syntax looks like for the most basic operations.

To get a more detailed introduction on Docker containers, see this other workshop on [Container workflows with Docker](https://pawseysc.github.io/container-workflows/).


### Downloading images and running containers

Let's cd into the relevant demo directory:

```
$ cd $TUTO/demos/lolcow_docker
```
{: .bash}

Let's download a Ubuntu container image, using `docker pull`:

```
$ sudo docker pull ubuntu:18.04
```
{: .bash}

```
18.04: Pulling from library/ubuntu
7ddbc47eeb70: Pull complete 
c1bbdc448b72: Pull complete 
8c3b70e39044: Pull complete 
45d437916d57: Pull complete 
Digest: sha256:6e9f67fa63b0323e9a1e587fd71c561ba48a034504fb804fd26fd8800039835d
Status: Downloaded newer image for ubuntu:18.04

```
{: .output}

Now, let's use this image via `docker run`:

```
$ sudo docker run ubuntu:18.04 cat /etc/os-release
```
{: .bash}

```
NAME="Ubuntu"
VERSION="18.04.3 LTS (Bionic Beaver)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 18.04.3 LTS"
VERSION_ID="18.04"
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
VERSION_CODENAME=bionic
UBUNTU_CODENAME=bionic
```
{: .output}

Similar to Singularity, the `pull` step could be skipped, as trying to run an image which is not locally available would trigger the download anyway.

To open up a shell in the image, use `docker run -it`:

```
$ sudo docker run -it ubuntu:18.04 bash
```
{: .bash}

```
root@dd1ca993f4ad:/#
```
{: .output}

Then type `exit`, or hit `Ctrl-D`, to leave the interactive shell.


As we mentioned above, lots of Docker defaults are about privileged runtime and container isolation.  Some extra syntax is required in order to achieve a container execution comparable to Singularity, *i.e.* with
* visibility of the host current working directory
* container working directory same as host one
* right user file ownership
* ability to pipe commands in the container

Long story short, this is what it takes:

```
$ sudo docker run --rm -v $(pwd):/data -w /data -u $(id -u):$(id -g) -i ubuntu:18.04 echo "Good Morning" >hello1.txt
$ ls -l hello1.txt
```
{: .bash}

```
-rw-r----- 1 ubuntu ubuntu 13 Nov  1 08:29 hello1.txt
```
{: .output}

Let's comment on the flags:
* `-v` is to bind mount host directories in the container
* `-w` is to set the container working directory
* `-u` is to set user/group in the container
* `-i` is to keep *STDIN* open in the container

What about the `--rm` flag? To respond to this, let's move on.


### Managing containers and images

By default, when containers exit, they remain cached in the system for potential future restart.  Have a look at a list of running and stopped containers with `docker ps -a` (remove `-a` to only list running ones):

```
$ sudo docker ps -a
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND                 CREATED             STATUS                       PORTS               NAMES
375a021f8674        ubuntu:18.04        "bash"                  52 seconds ago      Exited (0) 4 seconds ago                         reverent_volhard
6000f459c132        ubuntu:18.04        "cat /etc/os-release"   57 seconds ago      Exited (0) 55 seconds ago                        hungry_bhabha

```
{: .output}

It's possible to clean up cached, exited containers by means of `docker rm`; there's also an idiomatic way to clean all of them at once:

```
$ sudo docker rm $(sudo docker ps -qa)
```
{: .bash}

```
375a021f8674
6000f459c132
```
{: .output}

If I know in advance I won't need to re-run a container after it exits, I can use the runtime flag `--rm`, as in `docker run --rm`, to clean it up automatically, as we did in the example above.


Docker stores container images in a hidden directory under its own control.  To get the list of downloaded images use `docker images`:

```
$ sudo docker images
```
{: .bash}

```
REPOSITORY                        TAG                      IMAGE ID            CREATED             SIZE
ubuntu                            18.04                    775349758637        10 hours ago        64.2MB
```
{: .output}

If you don't need an image any more and want to clear up disk space, use `docker rmi` to remove it:

```
$ sudo docker rmi ubuntu:18.04
```
{: .bash}

```
Untagged: ubuntu:18.04
Untagged: ubuntu@sha256:6e9f67fa63b0323e9a1e587fd71c561ba48a034504fb804fd26fd8800039835d
Deleted: sha256:775349758637aff77bf85e2ff0597e86e3e859183ef0baba8b3e8fc8d3cba51c
Deleted: sha256:4fc26b0b0c6903db3b4fe96856034a1bd9411ed963a96c1bc8f03f18ee92ac2a
Deleted: sha256:b53837dafdd21f67e607ae642ce49d326b0c30b39734b6710c682a50a9f932bf
Deleted: sha256:565879c6effe6a013e0b2e492f182b40049f1c083fc582ef61e49a98dca23f7e
Deleted: sha256:cc967c529ced563b7746b663d98248bc571afdb3c012019d7f54d6c092793b8b
```
{: .output}


### A Dockerfile recipe

**Note**: the following sections on building and sharing container images with Docker are the same than in the dedicated episode.

It can be interesting to have an idea of how to build images with Docker.  In fact, as we mentioned earlier on, the layered image format of Docker can sometimes help in reducing image development time.  In addition, Docker images are quite universally compatible, as they can be run by Singularity, too.

We're going to build a very similar image to the one we built with Singularity.  The `Dockerfile` recipe file in the demo directory looks like:

```
FROM ubuntu:18.04

LABEL maintainer="Pawsey Supercomputing Centre"

RUN apt-get -y update && \
  apt-get -y install fortune cowsay lolcat

ENV PATH=/usr/games:$PATH

VOLUME /data
WORKDIR /data

CMD fortune | cowsay | lolcat
```
{: .source}

The directory where the the Dockerfile is stored is the so called the Docker **build context**.  Docker will include files in this directory in the build process and in the final image.  As a by-product, this will make the build process longer and the image larger, so that we want to include only those strictly required for the build, even none if possible.

Let's comment on the Docker instructions that appear in this Dockerfile.

* `FROM`: compulsory, it provides the starting image we will use to build our customised one;
* `LABEL`: used to add metadata information to the image, *e.g.* the maintainer, optional;
* `RUN`: this is the most used instruction, that allows to run most shell commands during the build.  Multiple `RUN` instructions are often found in a single Dockerfile;
* `ENV`: set environment variables that will persist at runtime in the container; **DO NOT** use `RUN export <..>` to this end, as the variable will be lost after the `RUN` step is completed;
* `VOLUME`: creates a mount point ready to be used for mounting external (*e.g.* host) volumes; creates the corresponding directory if not existing;
* `WORKDIR`: changes directory to the specified path; the last current directory in the build will be the working directory in the running container.  
  **Note**: if you use instead `RUN cd <..>`, the changed directory will only persist within that `RUN` instruction, and then be lost in subsequent build steps;
* `CMD`: specifies the default command to be executed with the container, in case no other command is provided.

More information on the Dockerfile syntax can be found at the [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).


### Layers in a container image

Note how the `RUN` instruction above is used to execute a sequence of commands to update the list of available packages and install a set of Linux packages.  
We have concatenated all these commands in one using the `&&` linux operator, and then the `\` symbol to break them into multiple lines for readability.

We could have used one `RUN` instruction per command, so why concatenating instead?  
Well, each `RUN` creates a distinct **layer** in the final image, increasing its size.  It is a good practice to use as few layers, and thus `RUN` instructions, as possible, to keep the image size smaller.


### Building the container image

Once the Dockerfile is ready, let us build the image with `docker build`:

```
$ sudo docker build -t lolcow:1Nov19 .
```
{: .bash}

In the command above, `.` is the location of the build context (*i.e.* the directory for the Dockerfile).  
The `-t` flag is used to specify the image name (compulsory) and tag (optional).

Any lowercase alphanumeric string can be used as image name; here we've used `lolcow`.  The image tag (following the colon) can be optionally used to maintain a set of different image versions on Docker Hub, and is a key feature in enabling reproducibility of your computations through containers; here we've used `1Nov19`.

Adding the prefix `<Your Docker Hub account>/` to the image name is also optional and allows to push the built image to your Docker Hub registry (see below). 

The complete format for the image name looks like: `<Your Docker Hub account ^>/<Image name>:<Image tag ^>`. `^`These are optional.

This is the output of our build:

```
Sending build context to Docker daemon  2.048kB
Step 1/7 : FROM ubuntu:18.04
 ---> 775349758637
Step 2/7 : LABEL maintainer="Pawsey Supercomputing Centre"
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

First, let us create a second tag for the image, that includes your Docker Account.  To this end we'll use `docker tag`:

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


### Sharing the Docker image as a single file

If you don't want to use an online registry to share your images, Docker allows you to convert them to a compressed `tar.gz` archive, which you can then share as any other large file, *e.g.* using `scp`, `rsync`, and other file transfer tools.  
For instance, this can be useful when needing to transfer or share images including proprietary software, amongst collaborators that own the appropriate license.

Use `docker save` to create the archive:

```
$ docker save -o lolcow_1Nov19.tar.gz lolcow:1Nov19
```
{: .bash}

After the transfer, use `docker load` to extract the image in a format that is usable by Docker:

```
$ docker load -i lolcow_1Nov19.tar.gz
```
{: .bash}

```
Loaded image: lolcow:1Nov19
```
{: .output}

**Note**: you need Docker to extract the image from the compressed archive, Singularity can't do it.
