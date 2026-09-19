---
title: "Setup Docker on your computer"
teaching: 10
exercises: 5
questions:
- How can I install Docker and verify that it works on my computer?
objectives:
- Install Docker on your own computer, on Windows, macOS or Linux
- Verify the installation by running a test container
- Create and verify a Docker Hub account and test publishing an image
keypoints:
- We will use Docker to build container images
- Installing Docker may require *admin*/*sudo* privileges, depending on your operating system and computer configuration
- The `docker run hello-world` command verifies that Docker can obtain an image and run a container
- A verified Docker Hub account is required for the later image-publishing exercise
---

### Why do I need this?

We will be using **Docker** to build Linux container images on your own computer.  Later on, we'll switch to **Singularity**, which is used on Setonix.  Docker is the main workshop software that you need to install *before* the session.  On Windows, the recommended Docker Desktop setup also uses WSL 2.

> ## Administrator permissions may be required
>
> Installing Docker, enabling WSL 2, or configuring required system features may require administrator permissions, depending on your operating system and computer configuration.  If you use a managed computer, contact your system administrator or IT support before the workshop.  Alternatively, use a personal computer on which you are permitted to install software.
{: .callout}

Pick the section below that matches your operating system.  In all cases, the goal is the same: get the `docker` command working, and confirm it by running a small test container called `hello-world`.


### 1. Windows

We recommend *Docker Desktop* with the *WSL 2* backend.  WSL stands for *Windows Subsystem for Linux*.  WSL 2 provides a Linux environment using a Linux kernel inside a lightweight virtual machine managed automatically by Windows.  Docker Desktop uses this environment to run Linux containers on Windows.

#### Step 1: Install WSL 2

Open *PowerShell* **as Administrator**, and run:

```
wsl --install
```
{: .bash}

This enables WSL 2 on your computer, and on most systems will also install *Ubuntu* as the default Linux distribution in the same step.  Restart your computer if prompted.

After restarting, open *PowerShell* again and check that WSL is working:

```
wsl --version
```
{: .bash}

This should print your installed WSL version, with no errors.

#### Step 2: Install Ubuntu on WSL

Check which Linux distributions are already installed:

```
wsl --list --verbose
```
{: .bash}

If `Ubuntu` is listed, you're done with this step.  If it isn't (or the list is empty), install it explicitly:

```
wsl --install -d Ubuntu
```
{: .bash}

The first time Ubuntu starts, it will ask you to create a Unix username and password — pick anything you like, you won't need them for this workshop.

#### Step 3: Install Docker Desktop

Download and install [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/).
During installation, use the **WSL 2 based engine** if prompted.  After installation, verify that Docker Desktop is configured to use WSL 2.  This training uses **Linux container images**, so Docker Desktop must run in Linux container mode, not Windows container mode.

Once installed, start Docker Desktop (it needs to be running in the background for the `docker` command to work).

#### Step 4: Check Docker

For consistency, we'll run the `docker` commands from **PowerShell**. Open *PowerShell* and run:

```
docker --version
```
{: .bash}

This confirms that the Docker command-line client is installed.  It does not confirm that the Docker engine is running.  To test the complete installation, run:

```
docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.

> ## Does `docker` also work inside the Ubuntu/WSL terminal?
>
> Yes. Docker Desktop exposes the `docker` command inside your WSL distributions too, through a setting called *WSL Integration* (*Settings → Resources → WSL Integration* in Docker Desktop), which is on by default for your default distribution.
>
> Under the hood, Docker Desktop runs inside its own `docker-desktop` WSL distribution, isolated from your Ubuntu one the same way any two WSL distributions are isolated from each other; it only talks to Ubuntu because WSL Integration is enabled for it. See Docker's [WSL 2 security in Docker Desktop](https://docs.docker.com/desktop/features/wsl/) for the full explanation.
{: .callout}


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

After installation, start Docker Desktop from the *Applications* folder or run:

```
open -a Docker
```
{: .bash}

Wait until Docker Desktop reports that the engine is running.  Then open *Terminal* and run:

```
docker --version
```
{: .bash}

This confirms that the Docker command-line client is installed.  It does not confirm that the Docker engine is running.  To test the complete installation, run:

```
docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect

> ## Additional check for Apple silicon Macs
>
> Setonix uses the `linux/amd64` platform for this training.  Apple silicon Macs use the `arm64` architecture, so Docker Desktop must use emulation to run the images used in the training.
>
> Run:
>
> ```
> docker run --rm --platform linux/amd64 alpine uname -m
> ```
> {: .bash}
>
> The expected output is:
>
> ```
> x86_64
> ```
> {: .output}
>
> If you see `x86_64`, Docker Desktop can pull and run `linux/amd64` containers on your Apple silicon Mac.
{: .callout}


### 3. Linux

The instructions below use *Ubuntu* as an example.  If you're on another distribution, Docker provides [installation instructions for several distributions](https://docs.docker.com/engine/install/).  For Ubuntu specifically, we recommend following the [Install using the `apt` repository](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository) section of the [Docker Engine on Ubuntu](https://docs.docker.com/engine/install/ubuntu/) guide for your version.

After installation, open a terminal and run:

```
docker --version
```
{: .bash}

This confirms that the Docker command-line client is installed.  It does not confirm that the Docker engine is running.  To test the complete installation, run:

```
sudo docker run hello-world
```
{: .bash}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.

> ## Running Docker without `sudo`
>
> On many Linux installations, Docker commands require `sudo`.  It is possible to run Docker without `sudo` by adding your account to the `docker` group, but membership in that group grants root-equivalent access to the system.
> This is **not required** for this training, so you can continue using `sudo docker ...`.
{: .callout}


### Final check

> ## Confirm your installation works
>
> Run the following command (`sudo docker run hello-world` on Linux, `docker run hello-world` in PowerShell on Windows, or in Terminal on macOS).  The exact output may vary depending on your Docker version and computer architecture, but it should include a message beginning with `Hello from Docker!`.
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
> > ```
> > {: .output}
> >
> > If you see a message starting with `Hello from Docker!` like this one, your installation is ready for the workshop — you're all set!
> {: .solution}
{: .challenge}


### Create and test a Docker Hub account

Later in the training, you will publish a container image to Docker Hub so that it can be pulled from Setonix with Singularity.  Create and verify your Docker Hub account before the workshop.

The steps below are based on Docker's official [Create a Docker account](https://docs.docker.com/accounts/individual/create-account/), [Create a repository](https://docs.docker.com/docker-hub/repos/create/), and [Push images to a repository](https://docs.docker.com/docker-hub/repos/manage/hub-images/push/) instructions.  You can follow those instructions directly if the Docker Hub interface or account requirements have changed.

1. Open the [Docker Hub sign-up page](https://hub.docker.com/signup/).
2. Create a free account using an email address, or continue with a supported external account.
3. Choose your Docker ID carefully.  Your Docker ID is the username used in Docker Hub image names, and it cannot be changed after the account is created.
4. Complete the account verification process.  You will not be able to sign in until the account has been verified.
5. Sign in to [Docker Hub](https://hub.docker.com/) and record your Docker ID.

In a terminal, assign your Docker ID to a shell variable.  Replace `<docker-id>` with your Docker ID and do not include the angle brackets:

```
$ DOCKER_ID="<docker-id>"
```
{: .bash}

#### Create a test repository

While signed in to Docker Hub:

1. Open **My Hub → Repositories**.
2. Select **Create repository**.
3. Select your personal Docker ID as the namespace.
4. Enter `first-image` as the repository name.
5. Set the repository visibility to **Public**.
6. Select **Create**.

The resulting repository name will be:

```
docker.io/<docker-id>/first-image
```
{: .output}

#### Test publishing an image

The earlier `docker run hello-world` test downloaded the `hello-world` image and ran a container from it.  Reuse that small image to verify that you can authenticate, tag an image for your Docker Hub namespace, and push it to the test repository.

Authenticate from the Docker client:

```
$ docker login docker.io
```
{: .bash}

Follow the authentication instructions shown by Docker.  Do not enter a password or access token directly as part of the command because doing so may record it in your shell history.

Create a new tag for the local `hello-world` image.  The new tag includes your Docker ID and the repository name:

```
$ docker tag hello-world:latest "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .bash}

Push the tagged image to Docker Hub:

```
$ docker push "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .bash}

After the push completes, open the `first-image` repository in Docker Hub and confirm that the `latest` tag is present.

To verify that Docker can retrieve the published image, remove its registry-qualified local tag, pull it from Docker Hub, and run it:

```
$ docker image rm "docker.io/${DOCKER_ID}/first-image:latest"
$ docker pull "docker.io/${DOCKER_ID}/first-image:latest"
$ docker run --rm "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .bash}

The output should begin with:

```
Hello from Docker!
```
{: .output}

This test publishes an existing small image rather than building a new one.  Image building, meaningful version tags, and publishing the training application are covered later in the workshop.

Do not upload proprietary, confidential, export-controlled, licensed, or otherwise restricted software to a public repository.


### If you run into problems

You can also refer to the official Docker documentation:
* [Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
* [macOS](https://docs.docker.com/desktop/setup/install/mac-install/)
* [Linux](https://docs.docker.com/engine/install/)


### Optional: further reading

If you'd like some additional background, Docker offers an introductory, self-paced workshop: [Getting Started with Docker](https://docs.docker.com/get-started/).
This is entirely optional — you do not need to complete it before attending this training.
