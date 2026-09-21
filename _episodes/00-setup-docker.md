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

We will use **Docker** to build Linux container images on your own computer. Later, we will switch to **SingularityCE**, referred to as **Singularity** throughout this training, to run the images on Setonix. Docker is the main workshop software that you need to install *before* the session. On Windows, the recommended Docker Desktop configuration uses the WSL 2 backend and an Ubuntu WSL environment for the workshop commands.

> ## Administrator permissions may be required
>
> Installing Docker, enabling WSL 2, or configuring required system features may require administrator permissions, depending on your operating system and computer configuration.  If you use a managed computer, contact your system administrator or IT support before the workshop.  Alternatively, use a personal computer on which you are permitted to install software.
{: .callout}

Pick the section below that matches your operating system.  In all cases, the goal is the same: get the `docker` command working, and confirm it by running a small test container called `hello-world`.


### 1. Windows

We recommend *Docker Desktop* with the *WSL 2* backend. WSL stands for *Windows Subsystem for Linux*. WSL 2 runs a Linux kernel inside a lightweight virtual machine managed automatically by Windows. Linux distributions such as Ubuntu run as isolated environments within WSL 2, while Docker Desktop uses the same WSL 2 infrastructure to run Linux containers on Windows.

Windows participants will use the Ubuntu WSL terminal for the local Docker exercises in this training. This provides the Bash shell and Linux command-line tools used throughout the lessons. Docker Desktop runs the Docker Engine in its own `docker-desktop` WSL distribution, while WSL Integration makes the `docker` command available from Ubuntu.

#### Step 1: Install WSL 2 and Ubuntu

Open *PowerShell* **as Administrator**. In the Windows setup instructions, `PS>` represents the PowerShell prompt. This prompt is a visual indicator and must not be typed as part of the command.

Run:

```powershell
PS> wsl --install
```
{: .source}

This enables WSL 2 on your computer and, on most systems, also installs *Ubuntu* as the default WSL Linux distribution. Restart your computer if prompted.

After restarting, open *PowerShell* again and check that WSL is working:

```powershell
PS> wsl --version
```
{: .source}

This should print your installed WSL version without errors.

Check which Linux distributions are installed and whether they use WSL 2:

```powershell
PS> wsl --list --verbose
```
{: .source}

The output should include an Ubuntu distribution. Note its exact name from the `NAME` column and confirm that the `VERSION` column shows `2`. For example:

```text
  NAME            STATE           VERSION
* Ubuntu-24.04    Stopped         2
```
{: .output}

If no Ubuntu distribution is listed, install it:

```powershell
PS> wsl --install -d Ubuntu
```
{: .source}

Then run `wsl --list --verbose` again and note the exact distribution name.

If the installed Ubuntu distribution is listed with version `1`, convert it to WSL 2 by using its exact name listed. For example, in this case `Ubuntu-24.04`:

```powershell
PS> wsl --set-version Ubuntu-24.04 2
```
{: .source}

Remember to replace `Ubuntu-24.04` with the exact name reported on your computer by the `wsl --list` command indicated above.

The first time Ubuntu starts, it will ask you to create a Unix username and password. These credentials belong only to your Ubuntu WSL environment and do not need to match your Windows or Pawsey credentials.

#### Step 2: Install Docker Desktop

Download and install [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/).

During installation, use the **WSL 2 based engine** if prompted. After installation, confirm that Docker Desktop is configured to use WSL 2. This training uses **Linux container images**, so Docker Desktop must run in Linux container mode, not Windows container mode.

Start Docker Desktop and leave it running in the background while using Docker commands.

#### Step 3: Enable Docker in Ubuntu

In Docker Desktop, open *Settings → Resources → WSL Integration*. Enable integration for the Ubuntu distribution and select *Apply*.

Docker Desktop runs the Docker Engine in its own `docker-desktop` WSL distribution. WSL Integration gives the Ubuntu environment access to the Docker command-line client and Docker Desktop engine without installing Docker Engine directly inside Ubuntu. See Docker's [WSL documentation](https://docs.docker.com/desktop/features/wsl/) for more information.

#### Step 4: Check Docker from Ubuntu

Open *Ubuntu* from the Windows Start menu. From this point onwards, Windows participants should run the local workshop commands in the Ubuntu terminal rather than PowerShell.

In the commands below, `$` represents the Bash prompt and must not be typed. Check the Docker client:

```bash
$ docker --version
```
{: .source}

This confirms that the Docker command-line client is available in Ubuntu. It does not confirm that the Docker Engine is running. To test the complete installation, run:

```bash
$ docker run hello-world
```
{: .source}

This downloads, if needed, and runs a small test image. See the [Final check](#final-check) section below for the output you should expect.

> ## Where should I keep the training files?
>
> On Windows, clone and work with the training repository inside the Ubuntu WSL filesystem, for example under your Ubuntu home directory (`~`). Avoid placing the working repository under `/mnt/c/` unless access from Windows is specifically required.
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

```bash
$ open -a Docker
```
{: .source}

Wait until Docker Desktop reports that the engine is running.  Then open *Terminal* and run:

```bash
$ docker --version
```
{: .source}

This confirms that the Docker command-line client is installed.  It does not confirm that the Docker engine is running.  To test the complete installation, run:

```bash
$ docker run hello-world
```
{: .source}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect

> ## Additional check for Apple silicon Macs
>
> Setonix uses the `linux/amd64` platform for this training.  Apple silicon Macs use the `arm64` architecture, so Docker Desktop must use emulation to run the images used in the training.
>
> Run:
>
> ```bash
> $ docker run --rm --platform linux/amd64 alpine uname -m
> ```
> {: .source}
>
> The expected output is:
>
> ```text
> x86_64
> ```
> {: .output}
>
> If you see `x86_64`, Docker Desktop can pull and run `linux/amd64` containers on your Apple silicon Mac.
{: .callout}


### 3. Linux

The instructions below use *Ubuntu* as an example.  If you're on another distribution, Docker provides [installation instructions for several distributions](https://docs.docker.com/engine/install/).  For Ubuntu specifically, we recommend following the [Install using the `apt` repository](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository) section of the [Docker Engine on Ubuntu](https://docs.docker.com/engine/install/ubuntu/) guide for your version.

After installation, open a terminal and run:

```bash
$ docker --version
```
{: .source}

This confirms that the Docker command-line client is installed.  It does not confirm that the Docker engine is running.  To test the complete installation, run:

```bash
$ sudo docker run hello-world
```
{: .source}

This downloads (if needed) and runs a small test image.  See the [Final check](#final-check) section below for the output you should expect.

> ## Running Docker without `sudo`
>
> On many Linux installations, Docker commands require `sudo`.  It is possible to run Docker without `sudo` by adding your account to the `docker` group, but membership in that group grants root-equivalent access to the system.
> This is **not required** for this training, so you can continue using `sudo docker ...`.
{: .callout}


### Final check

> ## Confirm your installation works
>
> Run the following command (`sudo docker run hello-world` on Linux, or `docker run hello-world` in the Ubuntu WSL terminal on Windows and in Terminal on macOS). The exact output may vary depending on your Docker version and computer architecture, but it should include a message beginning with `Hello from Docker!`.
>
> ```bash
> $ docker run hello-world
> ```
> {: .source}
>
> > ## Expected output
> >
> > ```text
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

Assign your Docker ID to a shell variable. Replace `<docker-id>` with your Docker ID and do not include the angle brackets. Run this command in Terminal on macOS or Linux, or in the Ubuntu WSL terminal on Windows:

```bash
$ DOCKER_ID="<docker-id>"
```
{: .source}

#### Create a test repository

While signed in to Docker Hub:

1. Open **My Hub → Repositories**.
2. Select **Create repository**.
3. Select your personal Docker ID as the namespace.
4. Enter `first-image` as the repository name.
5. Set the repository visibility to **Public**.
6. Select **Create**.

The resulting repository name will be:

```text
docker.io/<docker-id>/first-image
```
{: .output}

#### Test publishing an image

The earlier `docker run hello-world` test downloaded the `hello-world` image and ran a container from it.  Reuse that small image to verify that you can authenticate, tag an image for your Docker Hub namespace, and push it to the test repository.

Authenticate from the Docker client:

```bash
$ docker login docker.io
```
{: .source}

Follow the authentication instructions shown by Docker.  Do not enter a password or access token directly as part of the command because doing so may record it in your shell history.

Create a new tag for the local `hello-world` image.  The new tag includes your Docker ID and the repository name:

```bash
$ docker tag hello-world:latest "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .source}

Push the tagged image to Docker Hub:

```bash
$ docker push "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .source}

The Docker commands above use shell syntax that works in Bash and Zsh, including Bash in the Ubuntu WSL terminal on Windows.

After the push completes, open the `first-image` repository in Docker Hub and confirm that the `latest` tag is present.

To verify that Docker can retrieve the published image reference from Docker Hub, remove its registry-qualified local tag, pull it from Docker Hub, and run it:

```bash
$ docker image rm "docker.io/${DOCKER_ID}/first-image:latest"
$ docker pull "docker.io/${DOCKER_ID}/first-image:latest"
$ docker run --rm "docker.io/${DOCKER_ID}/first-image:latest"
```
{: .source}

The output should begin with:

```text
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
