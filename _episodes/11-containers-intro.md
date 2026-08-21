---
title: "Introduction to containers"
teaching: 10
exercises: 0
questions:
- What are containers for?
- Who is using containers in HPC ecosystems?
objectives:
- 'Define the term: "container" in contrast to "virtual machine"'
- 'Define other terms like: image and registry'
- Discuss when you would benefit from using containers in your workflow
keypoints:
- Containers enable you to package up an application and its dependencies.
- By using containers, you can better enforce reproducibility, portability and share-ability of your computational workflows.
---


### What is a container?

A container is a way of running one or more applications in an isolated software environment. The applications, tools, libraries and configuration they need are packaged together in a container image.

Containers are used to distribute software with its dependencies, avoid conflicts between software environments and make workflows easier to reproduce and move between compatible systems. They can be used across personal computers, cloud platforms and HPC systems.

### Containers vs Virtual Machines

If you have already used a Virtual Machine (VM), either on your own computer (for example, using VirtualBox) or through a cloud provider such as Azure, you're already familiar with some of the concepts needed to understand containers.

<!-- ![Containers vs. VMs]({{ page.root }}/fig/container_vs_vm.png) -->
<img class="img-responsive center-block" src="{{ page.root }}/fig/container_vs_vm.png" alt="Architecture of virtual machines and containers" style="width: 60%;"/>

VMs and containers provide isolation in different ways. A VM behaves like a complete computer running inside another computer, with its own guest operating system and kernel. A container isolates processes but does not boot its own kernel. Instead, **applications inside containers use the host system's Linux kernel**. (More generally, containers share the host system's kernel rather than running their own.)

The container image supplies the user-space environment required by the application, including applications, tools, libraries and a filesystem.

By sharing the host kernel, containers are generally:

* lighter weight to run (less CPU and memory usage)
* faster to start
* smaller in size (thus easier to transfer and share)

Container images are also typically built as specialised software environments for a particular application or workflow. This specialisation is a usage convention rather than an architectural difference: an image can contain several applications and tools, but they usually serve a common purpose.

Because containerised applications use the host kernel and CPU architecture, they must be compatible with the host system. For example, Linux containers require a compatible Linux kernel, and an image built for an `x86_64` CPU does not normally run on an `arm64` system, or vice versa. (Cross-architecture execution may be possible through emulation, but it is generally not appropriate for HPC workloads.)

### Why use containers?

There are a number of reasons for using containers in your daily work:

* Easier software installation and dependency management
  - Users can often run software directly from a container image provided by developers or software vendors, without installing anything themselves.
  - "I can't get this software stack to install in the cluster" is often what drives people to containers.

* Cross-system portability
  - Run the same software environment on your laptop, in the cloud and on HPC systems.
  - Greatly reduces the "it works on my machine" problem.

* Data reproducibility and software preservation
  - Helps ensure analyses can be repeated using the same software environment.
  - Preserve working software environments for months or years.
  - Useful when revisiting an old project after host operating systems, libraries or compilers have changed.

* Simplified collaboration
  - Share a complete software environment with collaborators
  - Then, avoid sharing lengthy installation instructions and configuration steps.

* Consistent testing environment
  - Test software in an environment that closely matches where it will run.
  - Reduce surprises caused by differences between systems.

A few examples of how containers are being used at Pawsey include:

* Bioinformatics workflows
* Machine Learning
* Python apps in radio astronomy
* RStudio & Jupyter Notebook sessions
* Webservers
* OpenFoam simulations
* Cloud workflows (via Singularity or Docker)
* HPC workflows (via Singularity)

Here's an overview of what a typical workflow looks like:

<!-- ![Container Workflow]({{ page.root }}/fig/container_lifecycle.png) -->
<img src="{{ page.root }}/fig/container_lifecycle.png" alt="Container Workflow" width="716" height="298"/>

### Terminology

An **image** is a file (or set of files) that contains an application together with its software dependencies, libraries, tools, run-time environment and filesystem. Images can be copied, shared, uploaded and downloaded.

A **container** is a running instance of an image. In other words, it is a process that has been started from an image. Multiple containers can be launched from the same image, just as the same application can be run multiple times with different inputs or options.

In abstract, an image corresponds to a file, whereas a container corresponds to a process.

A **registry** is a service that stores and distributes container images. Registries can be public (for example, *Docker Hub* or *Quay.io*) or private. Users can download images from registries and, where permitted, upload their own images for others to use.

To create an image, we need a recipe describing the software to install and how the image should be assembled. This recipe is called a **Definition File** (or **def file**) in the *Singularity* and *Apptainer* ecosystem, and a **Dockerfile** in the *Docker* ecosystem.


### Container engines

A number of tools are available to create, distribute and run containerised applications. Some of these will be covered throughout this tutorial:

* **Docker**: the most widely used container platform and image ecosystem. Docker is commonly used on personal computers, cloud systems and CI/CD platforms to build and distribute container images. Although Docker itself is not typically used directly on shared HPC systems, Docker images are commonly used as the starting point for HPC container workflows. See the extensive [docker documentation](https://docs.docker.com/) for more information.

* **Singularity**: a container engine designed for HPC environments, allowing users to run containers without requiring elevated privileges. Many existing HPC container workflows and images continue to use the Singularity name and image format. See the extensive [singularity documentation](https://docs.sylabs.io/guides/latest/user-guide/) for more information.

* **Apptainer**: an open-source continuation of the Singularity project. It provides the same core functionality and is becoming the standard open-source container engine on many HPC systems. Throughout this tutorial, most commands and concepts apply equally to Singularity and Apptainer. See the [apptainer documentation](https://apptainer.org/docs/user/main/) for more information.

Other container engines (not covered here) include:

* **Shifter/Sarus**: container runtimes designed for HPC systems with support for Docker-compatible images.
* **Charliecloud**: a lightweight container solution designed for HPC environments.
* **Enroot**: a lightweight container runtime developed by NVIDIA, commonly used for GPU-focused workloads.
* **Podman**: a daemonless, rootless container engine that is increasingly used as an alternative to Docker.

### Get ready for the hands-on

Before we start, let us ensure we have the required files to run the tutorials.

If you haven't done it already, download the following Github repo.  Then `cd` into it, and save the current directory into a variable named `TUTO` for later use.

```bash
$ cd ~
$ git clone https://github.com/PawseySC/singularity-containers
$ cd singularity-containers
$ export TUTO=$(pwd)
```
{: .source}

> ## Want to save time later in the tutorial?
>
> > ## Read this
> > Open a second terminal in the machine where you're running the tutorial, then run the script `pull_big_images.sh` to start downloading a few images that you'll require later:
> >
> > ```
> > $ cd $TUTO/demos
> > $ nohup bash ./pull_big_images.sh &
> > ```
> > {: .source}
> >
> > **Alternatively**, if you are running at Pawsey, *e.g.* on Zeus, submit this other script with Slurm instead:
> >
> > ```
> > $ cd $TUTO/demos
> > $ sbatch ./sbatch_pull_big_images.sh
> > ```
> > {: .source}
> >
> > This pull process will take at least one hour. Meanwhile, you'll be able to keep on going with this episode in your main terminal window.
> >
> {: .solution}
{: .challenge}


> ## Are you running on a shared HPC system?
>
> If you're running this tutorial on a shared system (*e.g.* on Zeus or Magnus at Pawsey), you should use one of the compute nodes rather than the login node.  You can get this setup by using an interactive scheduler allocation, for instance on Zeus with Slurm:
>
> ```
> $ salloc -n 1 -t 4:00:00
> ```
> {: .source}
>
> ```
> salloc: Granted job allocation 3453895
> salloc: Waiting for resource configuration
> salloc: Nodes z052 are ready for job
> ```
> {: .output}
{: .callout}
