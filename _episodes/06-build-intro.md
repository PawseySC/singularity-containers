---
title: "Build your own container image"
teaching: 20
exercises: 20
questions:
objectives:
- Learn what is a definition file (def file) and its basic syntax
- Learn how to build a container image and share it with others
- Learn the pros&cons of building with Singularity vs Docker
keypoints:
- "`sudo singularity build` is used to build images"
- "Use the remote builder with the flag `-r`, if you need to build images from a machine where you don't have sudo rights"
- "You can share you Singularity Image File with others, as you would do with any other (big) file"
- "`singularity push` can be used to upload images to a web registry (Sylabs account required)"
---


### What is a def file?

A *definition file*, or *def file*, is a recipe to build a container image with Singularity. It is basically a collection of the standard shell commands you would use to build your software through prompt; in addition, it contains Singularity-specific header lines that handle the build process. We will discuss these below with an example. Although there is no mandatory naming convention for def files, they are often characterised by the suffix `.def`.


### *Sudo* privileges with Singularity

Singularity does not allow for privileges escalation.  
In other words, if you are a standard user and you run `singularity`, any command inside the container will be run with the privileges of the standard user, *i.e.* without admin powers. If you try and `sudo` from inside the container you will get an error.  
On the other hand, if your user can run with *sudo*, and if you then decide to run Singularity as `sudo singularity`, then you will run any command from inside the container with admin powers.  
This design is what makes Singularity safe to run on HPC: users without admin rights are unable to escalate their privileges from inside the containers.

However, when building a container image you might need to install software using commands that require admin rights, *e.g.* `apt get` in Ubuntu/Debian or `yum` in Centos. To achieve this, you need to run `sudo singularity build`, implying that you need to carry out your build in a machine where you DO have admin rights.


### Building a basic container

Singularity can build container images in different formats. Let's focus on the Singularity Image Format, *i.e.* the one typically adopted to ship production-ready containers.

Let us cd into the appropriate directory:

```
$ cd $SC19/demos/07_lolcow
```
{: .bash}

Then, here is the def file we're going to use, `lolcow.def`:

```
BootStrap: docker
From: ubuntu:18.04

%post
    apt-get -y update
    apt-get -y install fortune cowsay lolcat

%environment
    export LC_ALL=C
    export PATH=/usr/games:$PATH

%labels
    Author Pawsey Supercomputing Centre
    Version v0.0.1

%help
    This is a demo container used to illustrate a def file.

%runscript
    fortune | cowsay | lolcat
```
{: .bash}

Let us build the image and run it first, then we'll comment on the contents of the def file. To this end we're using `sudo singularity build`, followed by the filename we decide to attribute to the container image, and then by the filename of the def file to be used:

```
$ sudo singularity build lolcow.sif lowcow.def
```
{: .bash}

```
INFO:    Starting build...
[..]
INFO:    Running post scriptlet
[..]
INFO:    Adding help info
INFO:    Adding labels
INFO:    Adding environment to container
INFO:    Adding runscript
INFO:    Creating SIF file...
INFO:    Build complete: lolcow.sif
```
{: .output}

Now, let us try and use the container simply as an executable:

```
$ ./lolcow.sif
```
{: .bash}

You will get something similar to this, hopefully just more colourful:

```
 _______________________________________
/ Have a place for everything and keep  \
| the thing somewhere else; this is not |
| advice, it is merely custom.          |
|                                       |
\ -- Mark Twain                         /
 ---------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}

Great, we've just containerised a cow that cites novelists! How did we achieve this?

The first line is `BootStrap: docker`.  
This tells Singularity how the image has to be initialised. `docker` means that we are going to start with a base image from Docker Hub. Another common way to bootstrap is using `library`, which will grab an image from the Sylabs Cloud. The image is specified in the next line, in this case `From: ubuntu:18.04`.  
Note how we started from Ubuntu 18.04 in Docker Hub, not Sylabs Cloud, as the former version has got a bit of a richer, more useful configuration.

Next is a section that start with the header `%post`. This is basically a sequence of commands to be executed to install packages in the image, in essence the same commands you would use for installation in a Linux box. Here we are ensuring we have an up-to-date list of packages, and then we are installing three Linux utilities.

The section `%environment` sets up environment variables that need to be defined at runtime rather than at build time. `LC_ALL` is required to avoid some locale warnings, whereas the `PATH` needs to be updated to reflect the location of the three utilities that we installed in the `%post` section.

Another section that is often useful can be defined by the header `%files`, like in:

```
%files
    <src-file> <dst-file>
```
{: .bash}

This section is used to copy files from the host, *i.e.* <src-file>, inside the container in the destination <dst-file>.


### Documenting the container image

The `%labels` section is used to add metadata to the container image. These can be then inspected by using

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
$ singularity inspect --deffile lolcow.sif 
```
{: .bash}


> ## Use the newly created container
> 
> You can use this new container using the same Singularity syntax we introduced > earlier on in this tutorial.
> 
> For instance, how would you run the command `fortune` from inside this > container?
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec lolcow.sif fortune
> > ```
> > {: .bash}
> > 
> > ```
> > Whenever you find that you are on the side of the majority, it is time
> > to reform.
> > 		-- Mark Twain
> > ```
> > {: .output}
> {: .solution}
> 
> Or, how would you open an interactive shell to explore the container?
> 
> > ## Solution
> > 
> > ```
> > $ singularity shell lolcow.sif 
> > ```
> > {: .bash}
> > 
> > ```
> > Singularity lolcow.sif:/some/dir> 
> > ```
> > {: .output}
> > 
> > Remember to close this session with `exit` or `Ctrl-D`.
> {: .solution}
{: .challenge}


### Run a container as an application

There's one section of the def file we haven't commented on yet. `%runscript` allows you to define a default command for the image. This command can then be used if you run the container as an executable:

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
$ singularity run -B $SC19/_episodes lolcow.sif
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


### Share your container image

Now that you've built your container image, you might want to run it on other systems, or share it with collaborators. 

The simplest way to achieve this is to remember that a SIF image is just a file, so .. you can transfer it across systems using Linux command line utilities like `scp` or `rsync`, or even graphical applications such as `Filezilla`.  
Just remember that images can be quite large, typically ranging from tens of MBs up to several GBs. For instance the *lolcow* image we created is about 70 MB.

If you want to keep the images publicly available, one option is to host them on the [**Sylabs Cloud Library**](https://cloud.sylabs.io), which is currently free upon signup. You can skip this and just follow the demo, if you don't want to signup.  
Once you create an account, you'll need to click on your account name on the top right of the page, select `Access Tokens`, then create a token, and copy it to the clipboard.  
Then you can configure the machine you're using for building container images, so that you can also push them to the Cloud Library:

```
$ singularity remote login
```
{: .bash}

```
Generate an API Key at https://cloud.sylabs.io/auth/tokens, and paste here:
API Key: 
```
{: .output}

Now paste the token you had copied to the clipboard end press `Enter`:

```
INFO:    API Key Verified!
```
{: .output}

You are now ready to push your image to the Cloud Library, *e.g.* via `singularity push`:

```
$ singularity push -U lolcow.sif library://<YOUR-SYLABS-USERNAME>/default/lolcow:30oct19
```
{: .bash}

```
WARNING: Skipping container verifying
 67.07 MiB / 67.07 MiB [==================================================================================================================================] 100.00% 6.59 MiB/s 10s
```
{: .output}

Note the use of the flag `-U` to allow pushing unsigned containers (see further down).  
Also note once again the format for the registry: <user>/<user-collection>/<name>:<tag>.

Finally, you (or other peers) are now able to pull your image from the Cloud Library:

```
$ singularity pull -U library://<YOUR-SYLABS-USERNAME>/default/lolcow:30oct19
```
{: .bash}

```
INFO:    Downloading library image
 67.07 MiB / 67.07 MiB [===================================================================================================================================] 100.00% 8.10 MiB/s 8s
WARNING: Skipping container verification
INFO:    Download complete: lolcow_30oct19.sif
```
{: .output}


### Remote build

What if you need to build an image from a system where you don't have admin privileges, *i.e.* you can't run commands with *sudo*?

Singularity offers the option to run build remotely, using the **Sylabs Remote Builder**; once again you will need a Sylabs account and a token to use this feature. If this is the case, just use `singularity build -r` to proceed with the remote build. Once finished, the image will be downloaded so that it's ready to use:

```
$ singularity build -r lolcow_remote.sif lolcow.def 
```
{: .bash}

```
INFO:    Remote "default" added.
INFO:    Authenticating with remote: default
INFO:    API Key Verified!
INFO:    Remote "default" now in use.
INFO:    Starting build...
[..]
INFO:    Running post scriptlet
[..]
INFO:    Adding help info
INFO:    Adding labels
INFO:    Adding environment to container
INFO:    Adding runscript
INFO:    Creating SIF file...
INFO:    Build complete: /tmp/image-699539270
WARNING: Skipping container verifying
 67.07 MiB / 67.07 MiB  100.00% 14.18 MiB/s 4s
```
{: .output}

At the time of writing, when using the Remote Builder you won't be able to use the `%files` header in the def file, to copy host files into the image.


### Other build options

The def file specification has a number of other interesting features, to know more about them you can visit the [Sylabs docs on def files](https://sylabs.io/guides/3.3/user-guide/definition_files.html).

In one of the following episodes we'll see how to use `%startscript` to configure the behaviour of containers running in background.

If you are in a development phase, where you don't know yet what you will include in your final container image, you can start with a *sandbox* image. This is a special type of image designed for development purposes, consisting not of a single file, but instead of a directory. To create one, run something like:

```
$ sudo singularity build --sandbox playbox/ docker://ubuntu:18.04
```
{: .bash}

Then to open it and play, run:

```
$ sudo singularity shell --writable playbox/
```
{: .bash}

More information on sandbox images can be found at the [Sylabs docs on building images](https://sylabs.io/guides/3.3/user-guide/build_a_container.html#creating-writable-sandbox-directories).

One last notable feature is the ability to use PGP keys to sign and verify container images. In this way, users of 3rd party containers can double check that the image they're running is bit-by-bit equivalent to the one that the author originally built, largely reducing the possibility to run containers infected by malware. you can find more on this topic at the [Sylabs docs on signing and verifying containers](https://sylabs.io/guides/3.3/user-guide/signNverify.html).


### Singularity *vs* Docker builds

We'll discuss how to build an image with Docker on a later episode. For now, let's just point out some of the advantages when building with one or the other tool. This will hopefully inform on which tool is best suited for you, depending on your specific context.

#### Singularity
* Single file image, can be handled as any other file
* Ability to sign/verify images
* Unambiguous container usage modes, via distinct keywords: `exec`, `shell`, `run`, `instance` (see episode on interactive containers)

#### Docker
* Image format can be run by all existing container engines
* Layered image format allows caching, for reduced build time during development phase


### Useful base images

At the time of writing, [Docker Hub](https://hub.docker.com) is the most popular web registry for general purpose container images. Therefore all images mentioned below are hosted in this registry.

#### CUDA

[nvidia/cuda](https://hub.docker.com/r/nvidia/cuda) has images to build GPU enabled applications. There are different image types for different needs. Tags containing `runtime` are suitable for binary applications that are ready to run; if you need to compile GPU code, pick tags containing `devel` instead. Different OS flavours are available, too.

#### MPI

As you can see in the episode on MPI applications, when containerising this type of software the MPI libraries in the image need to be ABI compatible with the MPI libraries in the host. The Pawsey Supercomputing Centre maintains some **MPICH** base images at [pawsey/mpi-base](https://hub.docker.com/r/pawsey/mpi-base), for building images that will run on our HPC systems.

#### Python

[python](https://hub.docker.com/_/python) hosts the official Python images. Different versions are available for some OS flavours. At the time of writing the default image tag corresponds to Python 3.8 on Debian 10. Smaller base images have tags ending with `-slim`.

[continuumio/miniconda3](https://hub.docker.com/r/continuumio/miniconda3) are images provided by the maintainers of the [Anaconda](https://anaconda.org) project. They ship with Python 3, as well as `pip` and `conda` to install and manage packages. At the time of writing, the most recent version is `4.7.12`, based on Python `3.7.4`.

If you need interactive Jupyter Notebooks, [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/) maintain a series of dedicated container images. Among others, there is the base SciPy image [jupyter/scipy-notebook](https://hub.docker.com/r/jupyter/scipy-notebook), the data science image [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook), and the machine learning image [jupyter/tensorflow-notebook](https://hub.docker.com/r/jupyter/tensorflow-notebook).

#### R

The [Rocker Project](https://www.rocker-project.org) maintains a number of good R base images. Of particular relevance is [rocker/tidyverse](https://hub.docker.com/r/rocker/tidyverse), which embeds the basic R distribution, an RStudio web-server installation and the [tidyverse](https://www.tidyverse.org) collection of packages for data science. At the time of writing, the most recent version is `3.6.1`.

Other more basic images are [rocker/r-ver](https://hub.docker.com/r/rocker/r-ver) (R only) and [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio) (R + RStudio).
