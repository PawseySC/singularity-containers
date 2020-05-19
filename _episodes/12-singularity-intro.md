---
title: "Basics of Singularity"
teaching: 20
exercises: 20
questions:
objectives:
- Learn how to download and run images
keypoints:
- Singularity can run both Singularity and Docker container images
- Execute commands in containers with `singularity exec`
- Open a shell in a container with `singularity shell`
- Download a container image in a selected location with `singularity pull`
---


### Get ready for the hands-on

Before we start, let us ensure we have got the required files to run the tutorials.

If you haven't done it already, download the following Github repo.  Then `cd` into it, define a couple of handy variables (see below), and finally `cd` into `demos/singularity`:

```
$ cd ~
$ git clone https://github.com/PawseySC/singularity-containers
$ cd singularity-containers
$ export TUTO=$(pwd)
$ cd demos/singularity
```
{: .bash}


> ## Want to save time later in the tutorial?  Read this
>
> Open a second terminal in the machine where you're running the tutorial, then run the script `pull_big_images.sh` to start downloading a few images that you'll require later:
>
> ```
> $ export TUTO=~/singularity-containers
> $ cd $TUTO/demos
> $ nohup bash ./pull_big_images.sh &
> ```
> {: .bash}
>
> This will take at least one hour. Meanwhile, you'll be able to keep on going with this episode in your main terminal window.
>
> One more thing: if you're running this tutorial on a shared system (*e.g.* on Zeus or Magnus at Pawsey), you should use one of the compute nodes rather than the login node.  You can get this setup by using an interactive scheduler allocation, for instance on Zeus with Slurm:
>
> ```
> $ salloc -n 1 -t 4:00:00
> ```
> {: .bash}
>
> ```
> salloc: Granted job allocation 3453895
> salloc: Waiting for resource configuration
> salloc: Nodes z052 are ready for job
> ```
> {: .output}
{: .callout}


### Singularity: a container engine for HPC

[Singularity](https://sylabs.io/singularity/) is developed and maintained by [Sylabs](https://sylabs.io), and was designed from scratch as a container engine for HPC applications, which is clearly reflected in some of its main features:

* *unprivileged* runtime: Singularity containers do not require the user to hold root privileges to run (the Singularity executable needs to be installed and owned by *root*, though);

* *root* privileges are required to build container images: users can build images on their personal laptops or workstations, on the cloud, or via a Remote Build service;

* *integration*, rather than *isolation*, by default: same user as host, same shell variables inherited by host, current directory bind mounted, communication ports available; as a result, launching a container requires a much simpler syntax than Docker;

* ability to run MPI enabled containers using host libraries;

* native execution of GPU enabled containers.

This tutorial assumes Singularity version 3.0 or higher. Version **3.3.0 or higher** is recommended as it offers a smoother, more bug-free experience.


### Container image formats

One of the differences between Docker and Singularity is the adopted format to store container images.

Docker adopts a layered format compliant with the *Open Containers Initiative* (OCI).  Each build command in the recipe file results in the creation of a distinct image layer.  These layers are cached during the build process, making them quite useful for development.  In fact, repeated build attempts that make use of the same layers will exploit the cache, thus reducing the overall build time.  On the other hand, shipping a container image is not straightforward, and requires either relying on a public registry, or compressing the image in a *tar* archive.

Since version 3.0, Singularity has developed the *Singularity Image Format* (SIF), a single file layout for container images, with extension `.sif`.  Among the benefits, an image is simply a very large file, and thus can be transferred and shipped as any other file.  Building on this single file format, a number of features have been developed, such as image signing and verification, and (more recently) image encryption.  A drawback of this approach is that during build time a progressive, incremental approach is not possible.

Interestingly, Singularity is able to download and run both types of images.

Note that Singularity versions prior to 3.0 used different image formats, characterised by the extensions `.simg` or `.sqsh`.  You can still find these around in the web; newer Singularity versions are still able to run them.


### Executing a simple command in a Singularity container

Running a command is done by means of `singularity exec`:

```
$ singularity exec library://ubuntu:18.04 cat /etc/os-release
```
{: .bash}

```
INFO:    Downloading library image

NAME="Ubuntu"
VERSION="18.04 LTS (Bionic Beaver)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 18.04 LTS"
VERSION_ID="18.04"
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
VERSION_CODENAME=bionic
UBUNTU_CODENAME=bionic
```
{: .output}

Here is what Singularity has just done:

* downloaded a Ubuntu image from the Cloud Library (this would be skipped if the image had been downloaded previously);
* stored it into the default cache directory;
* instantiated a container from that image;
* executed the command `cat /etc/os-release`.

Container images have a **name** and a **tag**, in this case `ubuntu` and `18.04`.  The tag can be omitted, in which case Singularity will default to a tag named `latest`.


> ## Using the *latest* tag
>
> The practice of using the `latest` tag can be handy for quick typing, but is dangerous when it comes to reproducibility of your workflow, as the *latest* image could change over time under the hood in the registry.
{: .callout}


Here Singularity pulled the image from an online image registry, as represented in this example by the prefix `library://`, that corresponds to the [**Sylabs Cloud Library**](https://cloud.sylabs.io).  Images in there are organised as: `<user>/<user collection>/<name>:<tag>`.  
In the example above we didn't we specify the **user**, `library`, and the **user collection**, `default`.  Why?  Because the specific case of `library/default/` can be omitted.  The full specification is used in the next example:

```
$ singularity exec library://library/default/ubuntu:18.04 echo "Hello World"
```
{: .bash}

```
Hello World
```
{: .output}

Here we are also experiencing image caching in action: the output has no more mention of the image being downloaded.


### Executing a command in a Docker container

Now let's try and download a Ubuntu container from the [**Docker Hub**](https://hub.docker.com), *i.e.* the main registry for Docker containers:

```
$ singularity exec docker://ubuntu:18.04 cat /etc/os-release
```
{: .bash}

```
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
Copying blob sha256:22e816666fd6516bccd19765947232debc14a5baf2418b2202fd67b3807b6b91
 25.45 MiB / 25.45 MiB [====================================================] 1s
Copying blob sha256:079b6d2a1e53c648abc48222c63809de745146c2ee8322a1b9e93703318290d6
 34.54 KiB / 34.54 KiB [====================================================] 0s
Copying blob sha256:11048ebae90883c19c9b20f003d5dd2f5bbf5b48556dabf06c8ea5c871c8debe
 849 B / 849 B [============================================================] 0s
Copying blob sha256:c58094023a2e61ef9388e283026c5d6a4b6ff6d10d4f626e866d38f061e79bb9
 162 B / 162 B [============================================================] 0s
Copying config sha256:6cd71496ca4e0cb2f834ca21c9b2110b258e9cdf09be47b54172ebbcf8232d3d
 2.42 KiB / 2.42 KiB [======================================================] 0s
Writing manifest to image destination
Storing signatures
INFO:    Creating SIF file...
INFO:    Build complete: /data/singularity/.singularity/cache/oci-tmp/a7b8b7b33e44b123d7f997bd4d3d0a59fafc63e203d17efedf09ff3f6f516152/ubuntu_18.04.sif

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

Rather than just downloading a SIF file, now there's more work for Singularity, as it has to both:

* download the various layers making up the image, and
* assemble them into a single SIF image file.

Note that, to point Singularity to Docker Hub, the prefix `docker://` is required.

Docker Hub organises images only by users (also called *repositories*), not by collections: `<reposityory>/<name>:<tag>`.  In the case of the Ubuntu image, the repository was `library` and could be omitted.


> ## What is the *latest* Ubuntu image from Docker Hub?
>
> Write down a Singularity command that prints the OS version through the *latest* Ubuntu image from Docker Hub.
>
> > ## Solution
> >
> > ```
> > $ singularity exec docker://ubuntu cat /etc/os-release
> > ```
> > {: .bash}
> >
> > ```
> > [..]
> > NAME="Ubuntu"
> > VERSION="20.04 LTS (Focal Fossa)"
> > [..]
> > ```
> > {: .output}
> >
> > It's version 20.04.
> {: .solution}
{: .challenge}


### Open up an interactive shell

Sometimes it can be useful to open a shell inside a container, rather than to execute commands, *e.g.* to inspect its contents.

Achieve this by using `singularity shell`:

```
$ singularity shell docker://ubuntu:18.04
```
{: .bash}

```
Singularity ubuntu_18.04.sif:/home/ubuntu/singularity-containers/demos/singularity>
```
{: .output}

Remember to type `exit`, or hit `Ctrl-D`, when you're done!


### Download and use images via SIF file names

All examples so far have identified container images using their registry name specification, *e.g.* `docker://ubuntu:18.04` or similar.

An alternative option to handle images is to download them to a known location, and then refer to their full directory path and file name.

Let's use `singularity pull` to save the image to a specified path (output might differ depending on the Singularity version you use):

```
$ singularity pull docker://ubuntu:18.04
```
{: .bash}

By default, the image is saved in the current directory:

```
$ ls
```
{: .bash}

```
ubuntu_18.04.sif
```
{: .output}

Then you can use this image file by:

```
$ singularity exec ./ubuntu_18.04.sif echo "Hello World"
```
{: .bash}

```
Hello World
```
{: .output}

You can specify the storage location with:

```
$ mkdir -p sif_lib
$ singularity pull --dir ~/path/to/sif/lib docker://library/ubuntu:18.04
```
{: .bash}

Being able to specify download locations allows you to keep the local set of images organised and tidy, by making use of a directory tree.  It also allows for easy sharing of images within your team in a shared resource.  In general, you will need to specify the location of the image upon execution:

```
$ singularity exec ~/path/to/sif/lib/ubuntu_18.04.sif echo "Hello Again"
```
{: .bash}

```
Hello Again
```
{: .output}


### Manage the image cache

When pulling images, Singularity stores images and blobs in a cache directory.

The default directory location for the image cache is `$HOME/.singularity/cache`.  This location can be inconvenient in shared resources such as HPC centres, where often the disk quota for the home directory is limited.  You can redefine the path to the cache dir by setting the variable `SINGULARITY_CACHEDIR`.

If you are running out of disk space, you can inspect the cache with this command (omit `-v` before Singularity version 3.4):

```
$ singularity cache list -v
```
{: .bash}

```
NAME                     DATE CREATED           SIZE             TYPE
ubuntu_18.04.sif         2020-05-14 08:30:05    54.62 MB         library
ubuntu_18.04.sif         2020-05-14 08:31:46    25.89 MB         oci
ubuntu_latest.sif        2020-05-14 08:32:35    27.77 MB         oci
23884877105a7ff84a9108   2020-05-14 08:31:36    26.69 MB         blob
2910811b6c4227c2f42aae   2020-05-14 08:31:38    0.85 kB          blob
357debb1f5fbc903e8d1e2   2020-05-14 08:31:42    0.81 kB          blob
36505266dcc64eeb1010bd   2020-05-14 08:31:40    0.16 kB          blob
6154df8ff9882934dc5bf2   2020-05-14 08:32:26    0.85 kB          blob
95c3f3755f37380edb2f8f   2020-05-14 08:32:29    2.48 kB          blob
96878229af8adf91bcbf11   2020-05-14 08:32:29    0.81 kB          blob
a0a91ba64d5bae825e23f0   2020-05-14 08:31:41    2.48 kB          blob
bc38caa0f5b94141276220   2020-05-14 08:31:37    35.37 kB         blob
d51af753c3d3a984351448   2020-05-14 08:32:24    28.56 MB         blob
fc878cd0a91c7bece56f66   2020-05-14 08:32:25    32.30 kB         blob
fee5db0ff82f7aa5ace634   2020-05-14 08:32:27    0.16 kB          blob

There are 3 container file(s) using 108.28 MB and 12 oci blob file(s) using 55.32 MB of space
Total space used: 163.60 MB
```
{: .output}

we are not going to clean the cache in this tutorial, as cached images will turn out useful later on.  Let us just perform a dry-run using the `-n` option:

```
$ singularity cache clean -n
```
{: .bash}

```
User requested a dry run. Not actually deleting any data!
Removing /home/ubuntu/.singularity/cache/library
Removing /home/ubuntu/.singularity/cache/oci-tmp
Removing /home/ubuntu/.singularity/cache/shub
Removing /home/ubuntu/.singularity/cache/oci
Removing /home/ubuntu/.singularity/cache/net
Removing /home/ubuntu/.singularity/cache/oras
```
{: .output}

If we really wanted to wipe the cache, we would need to use the `-f` flag instead (or, before Singularity version 3.4, the `-a` flag).


> ## Contextual help on Singularity commands
>
> Use `singularity help`, optionally followed by a command name, to print help information on features and options.
{: .callout}


### Popular registries (*aka* image libraries)

At the time of writing, Docker Hub hosts a much wider selection of container images than Sylabs Cloud.  This includes Linux distributions, Python and R deployments, as well as a big variety of applications.

Bioinformaticians should keep in mind another container registry, [Quay](https://quay.io) by Red Hat, that hosts thousands of applications in this domain of science.  These mostly come out of the [Biocontainers](https://biocontainers.pro) project, that aims to provide automated container builds of all of the packages made available through [Bioconda](https://bioconda.github.io).

Nvidia maintains the [Nvidia GPU Cloud (NGC)](https://ngc.nvidia.com), hosting an increasing number of containerised applications optimised to run on GPUs.

Right now, the Sylabs Cloud Library does not contain a large number of images.  Still, it can turn useful for storing container images requiring features that are specific to Singularity (we will see some in the next episodes).


> ## Pull and run a Python container ##
>
> How would you pull the following container image from Docker Hub, `python:3-slim`?
>
> Once you've pulled it, enquire the Python version inside the container by running `python --version`.
>
> > ## Solution
> >
> > Pull:
> >
> > ```
> > $ singularity pull docker://python:3-slim
> > ```
> > {: .bash}
> >
> > Get Python version:
> >
> > ```
> > $ singularity exec ./python_3-slim.sif python --version
> > ```
> > {: .bash}
> >
> > ```
> > Python 3.8.0
> > ```
> > {: .output}
> {: .solution}
{: .challenge}
