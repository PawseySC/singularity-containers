---
title: "Introduction to Containers"
teaching: 10
exercises: 0
questions:
objectives:
- Learn what are containers and how they are beneficial to HPC & Cloud computing
---


### Containers vs Virtual Machines

Many of you have probably used a VM, so you're actually already familiar with some of the concepts of a container.

![Containers vs. VMs]({{ page.root }}/fig/container_vs_vm.png)

The key difference here is that VMs virtualise **hardware** while containers virtualise **operating systems**.  There are other differences (and benefits)

* Containers are lighter weight (less CPU and memory usage, faster start-up times)
* More portable
* Modular (can easily combine multiple containers that work together)


### Terminology

An **image** is a file (or set of files) that contains the application and all its dependencies, libraries, run-time systems, etc. required to run. You can copy images around, upload them, download them etc.

A **container** is an instantiation of an image.  That is, it's a running process spawned out of an image.  You can run multiple containers from the same image, much like you might run the same application with different options or arguments.

In general, an image corresponds to a file, a container corresponds to a process.

A **registry** is a server application where images are stored and can be accessed by users. It can be public (e.g. Docker Hub) or private.

To build an image we need a recipe. A recipe file is called a **Definition File**, or **def file**, in the *Singularity* jargon and a **Dockerfile** in the *Docker* world.


### Containers and your workflow

There are a number of reasons for using containers in your daily work:

* Data reproducibility/provenance
* Cross-platform portability
* Simplified collaboration
* Simplified software dependencies and management
* Scalability
* Consistent testing environment

A few examples of how containers are being used at Pawsey include:

* Bioinformatics workflows
* Machine Learning 
* Python apps in radio astronomy
* RStudio & Jupyter Notebook sessions
* Webservers
* OpenFoam simulations
* Cloud workflows (via Singularity or Docker)
* HPC workflows (via Singularity or Shifter)

Here's an overview of what a typical workflow looks like:

![Docker Workflow]({{ page.root }}/fig/container_lifecycle.png)


### Container engines

A number of tools are available to create, deploy and run containerised applications. Some of these will be covered throughout this tutorial:

* **Docker**: the first engine to gain popularity, still widely used in the IT industry. Not very suitable for HPC as it requires *root* privileges to run. We'll use it mostly to build container images.

* **Singularity**: a simple, powerful container engine for the HPC world. The main focus of this workshop.

* **Shifter/Sarus**: a Docker-compatible container engine, suitable for HPC. Can run containers, cannot build them.

* **Podman**: a root-less alternative to Docker. Misses on some features for HPC.
