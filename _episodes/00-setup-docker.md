---
title: "Setup Docker on your computer"
teaching: 10
exercises: 5
questions:
objectives:
- Install Docker on your own computer, on Windows, macOS or Linux
- Verify the installation by running a test container
keypoints:
- We will use Docker to build container images
- Installing Docker requires *admin*/*sudo* privileges on your machine
- The 'docker run hello-world' command is the quickest way to check that your installation works
---

### Why do I need this?

We will be using **Docker** to build container images on your own computer.  Later on, we'll switch to **Singularity/Apptainer**, which is used on Setonix.  So, Docker is the one piece of software you need to install *before* the session.

> ## Administrator permissions required
>
> Installing WSL and Docker requires administrator permissions on your computer.  If you do not have these permissions, contact your system administrator or IT support before the workshop.  Alternatively, use a personal laptop on which you can install software.
{: .callout}

Pick the section below that matches your operating system.  In all cases, the goal is the same: get the `docker` command working, and confirm it by running a small test container called `hello-world`.


### 1. Windows

We recommend *Docker Desktop* with the *WSL 2* backend (WSL stands for *Windows Subsystem for Linux*, and lets you run a real Linux environment alongside Windows).

#### Step 1: Install WSL 2

Open *PowerShell* **as Administrator**, and run:

```
wsl --install
```
{: .bash}

This installs WSL and a default Linux distribution.  Restart your computer if prompted.

After restarting, open *PowerShell* again and check that WSL is working:

```
wsl --version
```
{: .bash}

This should print your installed WSL version, with no errors.

#### Step 2: Install Docker Desktop

Download and install [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/).  
During installation, make sure you select the **WSL 2 based engine** when prompted.

Once installed, start Docker Desktop (it needs to be running in the background for the `docker` command to work).

#### Step 3: Check Docker

Open *PowerShell* (or your terminal of choice) and run:

```
docker --version
```
{: .bash}

This checks that Docker is installed and that the command is available. Then run the test container:

```
docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.


### 2. macOS

Install [Docker Desktop for Mac](https://docs.docker.com/desktop/setup/install/mac-install/).

> ## Pick the right version
>
> Make sure you download the version that matches your Mac's processor:
> * **Apple silicon** (M1, M2, M3, M4, ...)
> * **Intel**
>
> If you're not sure which one you have, check via the Apple menu → *About This Mac*.
{: .callout}

After installation, start Docker Desktop, then open *Terminal* and run:

```
docker --version
```
{: .bash}

This checks that Docker is installed and that the command is available. Then run the test container:

```
docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.


### 3. Linux

The instructions below use *Ubuntu* as an example.  If you're on another distribution, Docker provides [installation instructions for several distributions](https://docs.docker.com/engine/install/).  For Ubuntu specifically, we recommend following the [Install using the `apt` repository](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository) section of the [Docker Engine on Ubuntu](https://docs.docker.com/engine/install/ubuntu/) guide for your version.

After installation, open a terminal and run:

```
docker --version
```
{: .bash}

This checks that Docker is installed and that the command is available. Then run the test container:

```
sudo docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.

> ## Running Docker without `sudo`
>
> Depending on your setup, you can configure Docker so that you don't need `sudo` for every command (by adding your user to the `docker` group).  
> This is **not required** for this training — `sudo docker ...` works just as well — so feel free to skip it for now.
{: .callout}


### Final check

> ## Confirm your installation works
>
> Run the following command (`sudo docker run hello-world` on Linux, `docker run hello-world` on Windows/macOS), and check that your output looks like the one below.
>
> ```
> $ docker run hello-world
> ```
> {: .bash}
>
> > ## Expected output
> >
> > ```
> > Hello from Docker!
> > This message shows that your installation appears to be working correctly.
> >
> > To generate this message, Docker took the following steps:
> >  1. The Docker client contacted the Docker daemon.
> >  2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
> >     (amd64)
> >  3. The Docker daemon created a new container from that image which runs the
> >     executable that produces the output you are currently reading.
> >  4. The Docker daemon streamed that output to the Docker client, which sent it
> >     to your terminal.
> >
> > To try something more ambitious, you can run an Ubuntu container with:
> >  $ docker run -it ubuntu bash
> >
> > Share images, automate workflows, and more with a free Docker ID:
> >  https://hub.docker.com/
> >
> > For more examples and ideas, visit:
> >  https://docs.docker.com/get-started/
> > ```
> > {: .output}
> >
> > If you see a message starting with `Hello from Docker!` like this one, your installation is ready for the workshop — you're all set!
> {: .solution}
{: .challenge}


### If you run into problems

You can also refer to the official Docker documentation:
* [Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
* [macOS](https://docs.docker.com/desktop/setup/install/mac-install/)
* [Linux](https://docs.docker.com/engine/install/)


### Optional: further reading

If you'd like some additional background, Docker offers an introductory, self-paced workshop: [Getting Started with Docker](https://docs.docker.com/get-started/).  
This is entirely optional — you do not need to complete it before attending this training.
