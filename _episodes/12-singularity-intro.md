---
title: "Basic Use of Containers (with Singularity)"
teaching: 35
exercises: 25
questions:
- "How can I find and run an existing container image with Singularity?"
- "How can I save and organise container images as SIF files?"
- "What is the difference between `singularity run`, `exec`, and `shell`?"
- "How can I execute pipelines and other shell expressions inside a container?"
objectives:
- "Identify the components of a Docker/OCI image reference"
- "Download and organise an image as a local SIF file"
- "Run predefined and user-selected commands from a container"
- "Inspect the software environment packaged inside a container"
- "Run pipelines and other shell expressions inside a container with `bash -c`"
- "Open an interactive shell inside a container"
keypoints:
- "Singularity can retrieve Docker/OCI images from compatible registries without using the Docker engine"
- "`singularity run` executes the action defined by the image publisher"
- "`singularity exec` executes a command selected by the user"
- "`singularity shell` opens an interactive shell for inspection and troubleshooting"
- "Use `bash -c` when commands executed with `singularity exec` require shell built-ins, pipelines, redirections, or other shell syntax"
- "Use `singularity pull` to create explicitly named SIF files that you can organise and reuse"
- "For images to be downloaded, prefer versioned image tags over `latest` when they are available"
---

### Get ready for the hands-on

Before we start, let us ensure we have got the required files to run the tutorials.

If you haven't done so already, move to a suitable working directory and download the following GitHub repository. On Pawsey systems, use your scratch directory; on other HPC or cloud systems, use the equivalent working directory recommended by the system administrators.

```bash
$ cd "$MYSCRATCH"    # On Pawsey systems
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
$ cd "$TUTO"
```
{: .source}

Now `cd` to the working directory. In this case:
```bash
$ cd demos/basic_use
$ pwd
```

The working directory should be something like:
```text
/path/to/scratch/singularity-containers/demos/basic_use
```
{: .source}



<div class="panel panel-warning">
  <div class="panel-heading">
    <strong>Content update required — start</strong><br>
    Update the following hands-on instructions
  </div>

  <div class="panel-body" markdown="1">

> ## Want to save time later in the tutorial?
>
> > ## Read this
> > Open a second terminal in the machine where you're running the tutorial, then run the script `pull_big_images.sh` to start downloading a few images that you'll require later:
> >
> > ```
> > $ cd $TUTO/demos
> > $ nohup bash ./pull_big_images.sh &
> > ```
> > {: .bash}
> >
> > **In alternative**, if you are running at Pawsey, *e.g.* on Zeus, submit this other script with Slurm instead:
> >
> > ```
> > $ cd $TUTO/demos
> > $ sbatch ./sbatch_pull_big_images.sh
> > ```
> > {: .bash}
> >
> > This pull process will take at least one hour. Meanwhile, you'll be able to keep on going with this episode in your main terminal window.
> >
> {: .solution}
{: .challenge}


> ## Are you running on a shared HPC system?
>
> If you're running this tutorial on a shared system (*e.g.* on Setonix at Pawsey), you should use one of the compute nodes rather than the login node. You can get this setup by using an interactive scheduler allocation, for instance on Setonix with Slurm:
>
> ```
> $ salloc -N 1 -n 1 -c 1 --reservation=ContainersTraining -t 4:00:00
> ```
> {: .bash}
>
> ```
> salloc: Granted job allocation 3453895
> salloc: Waiting for resource configuration
> salloc: Nodes nid002604 are ready for job
> ```
> {: .output}
{: .callout}

  </div>

  <div class="panel-footer">
    <strong>Content update required — end</strong>
  </div>
</div>


### Singularity: the container engine used in this training

Singularity is a container engine designed for shared HPC environments. It allows users to run containers without requiring elevated privileges and integrates containerised applications with host filesystems, schedulers, networks, and HPC hardware.

Throughout this training, we use Singularity to run containers. Many of the images used with Singularity were originally built and published through the Docker/OCI ecosystem. Singularity can retrieve these images directly from compatible registries, convert them to the Singularity Image Format (SIF), and run them without requiring the Docker engine.

Most of the commands and concepts in this episode also apply to Apptainer, an open-source continuation of the Singularity project.

### Exploring scientific images on Docker Hub

As discussed in the introductory episode, container registries store and distribute container images. Before building your own image, check whether the application developers, a software vendor, or another trusted organisation already provides a suitable image.

[Docker Hub](https://hub.docker.com/) is a widely used public container registry. It contains images published by software vendors, open-source projects, organisations, and individual users.

Before using our first image, let us explore how scientific software is presented on Docker Hub.

> ## Explore Docker Hub
>
> Open [Docker Hub](https://hub.docker.com/) in a web browser and search for each of the following image repositories:
>
> ```text
> rocker/rstudio
> tensorflow/tensorflow
> opencfd/openfoam-default
> ```
> {: .output}
>
> These repositories provide images for different types of scientific work:
>
> - `rocker/rstudio` provides an RStudio Server environment for statistical computing and data analysis.
> - `tensorflow/tensorflow` provides the TensorFlow machine-learning framework, with different image variants for CPU, GPU, development, and interactive environments.
> - `opencfd/openfoam-default` provides OpenFOAM applications, runtime libraries, source code, development tools, and tutorial cases for computational fluid dynamics.
>
> For each repository, inspect:
>
> - the organisation or user that published the image
> - the image description and documentation
> - the available tags
> - when the image was last updated
> - the supported CPU architectures
> - the approximate compressed image size
> - whether source files or build instructions are linked
> - whether the repository provides different image variants
>
> Do not download any of these images yet. Some scientific application images are large, and we are only exploring how images and their metadata are presented in Docker Hub.
>
> Were there other repositories with similar names in the search results? What information would help you decide which publisher to trust?
{: .challenge}

An image repository can provide many related images through its tags. A tag may identify:

- an application version
- a base operating-system version
- CPU or GPU support
- a minimal or extended software environment
- a stable release or a development build

Consequently, selecting a container image requires more than finding a familiar application name. You must also identify the publisher and select an appropriate image tag.

Docker Hub contains images from many publishers. Finding an image with the desired name does not automatically mean that it is trustworthy, maintained, compatible with your target system, or suitable for research use. When evaluating an image, consider:

- who published it
- whether it is maintained by the application developers or another trusted organisation
- when it was last updated
- which tags are available
- whether documentation or build instructions are provided
- whether the image is built for the operating system and CPU architecture of the target system
- for GPU-enabled applications, whether the image supports the GPU vendor and architecture, and is compatible with the GPU software environment on the target system

### Finding the image used in this episode

The scientific images that we have explored provide realistic examples, but some of them are large or require additional configuration. For our first commands, we will use a smaller and simpler teaching image.

Search Docker Hub for:

```text
sylabsio/lolcow
```
{: .output}

Open the repository named `sylabsio/lolcow`.

The repository name contains:

- `sylabsio`: the publisher's namespace
- `lolcow`: the image repository

The complete tagged image name is:

```text
sylabsio/lolcow:latest
```
{: .output}

The `latest` component is the image tag. We will discuss the meaning and limitations of this tag later in the episode.

> ## A teaching image
>
> We use `sylabsio/lolcow` because it produces an immediate and visually distinctive result without requiring input data or application-specific knowledge.
>
> It should not be interpreted as a recommendation for research or production workloads. Its role is to help us learn the basic Singularity commands before working with larger scientific application images.
{: .callout}

> ## Running this training on Setonix
>
> The following setup is specific to Setonix. If you are completing this training on another system, use the Singularity or Apptainer installation provided by that system.
>
> Before continuing, confirm that you are working inside an interactive `salloc` session on a Setonix compute node, as described in the hands-on preparation section above. **Do not run the following container commands directly on a login node.**
>
> You can confirm the current node with:
>
> ```bash
> $ hostname
> ```
> {: .source}
>
> A Setonix compute-node hostname begins with `nid`.
>
> Once you are on a compute node, load the Singularity module:
>
> ```bash
> $ module load singularity/4.1.0-nohost
> ```
> {: .source}
>
> On Setonix, multiple Singularity module variants are available. These variants configure different levels of integration between the container and the host software environment, including support for MPI, GPUs, and Slurm.
>
> The `nohost` suffix identifies the module variant intended for containers that should remain isolated from the host software environment. This avoids introducing host MPI libraries or other specialised host integrations that are not required by these introductory examples.
>
> Confirm that the expected Singularity version is available after loading the module:
>
> ```bash
> $ module list
> ```
> {: .source}
>
> The correct module should be listed as loaded:
>
> ```text
> Currently Loaded Modules:
  ...
  15) singularity/4.1.0-nohost
> ```
> {: .output}
>
> Check the version directly using the now available `singularity` command:
>
> ```bash
> $ singularity --version
> ```
> {: .source}
>
> The output should identify SingularityCE version 4.1.0:
>
> ```text
> singularity-ce version 4.1.0
> ```
> {: .output}
>
> You only need to load the module once in each terminal session. If you open another terminal, start a new login session, or submit a batch job, load the module again in that environment before invoking `singularity`.
>
> The available module versions and variants can change when the Setonix software environment is updated. To see the currently available modules, run:
>
> ```bash
> $ module avail singularity
> ```
> {: .source}
>
> For descriptions of the available module variants, see the
> [Singularity documentation in the Pawsey User Support Documentation](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925894/Singularity).
{: .callout}

### Downloading an image as a SIF file

Singularity uses the Singularity Image Format (SIF) for its native container images. A SIF image is a single file that can be copied, renamed, moved, and stored like any other file.

For images that you intend to keep and reuse, we recommend storing named SIF files in an organised personal or project image library.

First, create a directory for your Singularity images:

```bash
$ mkdir -p "$MYSOFTWARE/singularity/images"
```
{: .source}

The source image that we found on Docker Hub is identified by:

```text
docker://docker.io/sylabsio/lolcow:latest
```
{: .output}

The components of this image reference are:

- `docker://`: tells Singularity to retrieve a Docker/OCI image through a compatible container registry
- `docker.io`: the registry hostname that explicitly identifies Docker Hub
- `sylabsio`: the namespace of the organisation or user that published the image
- `lolcow`: the image repository name
- `latest`: the image tag

The `docker://` prefix does not instruct Singularity to start Docker. Singularity communicates directly with the registry and processes the Docker/OCI image manifest and filesystem layers. Docker does not need to be installed or running on the system.

First, create a directory (your local library) where to keep your singularity images:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```
{: .source}

Then, use `singularity pull` to download the source image, convert it to SIF, and save it with an explicit filename:

```bash
$ SINGULARITY_IMAGE="${MY_LOCAL_LIBRARY}/lolcow--latest.sif"
$ singularity pull "$SINGULARITY_IMAGE" docker://docker.io/sylabsio/lolcow:latest
```
{: .source}

The first argument after `pull` is the SIF file that Singularity creates (in this case the name that SINGULARITY_IMAGE has):

```text
${MYSOFTWARE}/singularity/images/lolcow--latest.sif
```
{: .output}

The second argument is the source image reference:

```text
docker://docker.io/sylabsio/lolcow:latest
```
{: .output}

Check that the SIF file was created:

```bash
$ ls -lh "$SINGULARITY_IMAGE"
```
{: .source}

A filesystem directory containing SIF files is not technically a container registry. It is an organised local collection, or image library, that you manage yourself.

<details markdown="1">
<summary class="alert alert-info"><strong>Optional: Docker Hub shorthand</strong></summary>

Docker Hub is the default registry for `docker://` references. Therefore:

```text
docker://sylabsio/lolcow:latest
```
{: .output}

is equivalent to:

```text
docker://docker.io/sylabsio/lolcow:latest
```
{: .output}

This tutorial uses the explicit `docker.io` hostname to make the registry visible. You will commonly encounter the abbreviated form in documentation and existing workflows.

</details>

> ## `pull` or `build`?
>
> Use `singularity pull` when obtaining an existing image from a registry:
>
> ```bash
> $ singularity pull output.sif docker://registry/namespace/image:tag
> ```
> {: .source}
>
> `singularity build` would have also worked for this example, but it is a more general command. It can also create a SIF image from a registry reference, but it is principally introduced later when building or customising images from definition files.
{: .callout}

### Running a container's default action

The image is now available as a SIF file in your personal image library. Define a variable containing its path so that it can be referenced more conveniently in subsequent commands:

```bash
$ SINGULARITY_IMAGE="${MYSOFTWARE}/singularity/images/lolcow--latest.sif"
```
{: .source}

Singularity provides three main commands for running or interacting with a container image:

- `singularity run` executes the image's predefined run action.
- `singularity exec` executes a command selected by the user.
- `singularity shell` opens an interactive shell inside the container environment.

We will use all three commands in this episode. First, use `singularity run` to execute the default action provided by the `lolcow` image:

```bash
$ singularity run "$SINGULARITY_IMAGE"
```
{: .source}

```text
 ______________________________
< Fri Sep 4 16:54:42 AWST 2026 >
 ------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}

The exact message and colours will vary because the image generates the output dynamically.

The image's default action combines several command-line programs:

- `date` prints the current date and time.
- `cowsay` displays that text in a speech bubble above an ASCII-art figure.
- `lolcat` adds changing terminal colours.

Conceptually, the image runs a pipeline similar to:

```bash
date | cowsay | lolcat
```
{: .source}

These programs were not developed specifically for Singularity, and they can be installed directly on many Linux distributions. In this example, however, they are provided by the container image and do not need to be installed on the host system. (BTW, `lolcow` is just the name of the image but not an existing command or script in the image or anywhere.)

You can check whether the commands are available directly on the host:

```bash
$ command -v date cowsay lolcat
```
{: .source}

If a command is not available, `command -v` does not print a path for it. Regardless of whether any of these programs happen to be installed on the host, the copies packaged inside the container are available when the container is used.

This illustrates an important benefit of containers: the required applications are supplied by the image instead of depending on software installed separately on each host system.

The `run` command starts a container from the SIF image and executes the image's predefined **runscript**. The runscript is configured by the image publisher when the image is built.

The general form of the command is:

```bash
$ singularity run IMAGE [ARGUMENTS...]
```
{: .source}

Running an image does not necessarily open an interactive session. The action performed by `singularity run` depends on the runscript defined in that particular image. For `lolcow`, the runscript generates a random message and displays it using an ASCII-art cow.

In the following sections, we will use:

- `singularity exec` to choose a particular command to run from the image
- `singularity shell` to explore the container environment interactively

<details markdown="1">
<summary class="alert alert-info"><strong>Optional: Running an image using an online registry reference</strong></summary>

Creating a named SIF file with `singularity pull` gives you control over where the image is stored and how it is named. This is the recommended approach for images that you intend to manage and reuse.

For a quick test, Singularity can also retrieve and run an image using its reference in an online registry:

```bash
$ singularity run docker://docker.io/sylabsio/lolcow:latest
```
{: .source}

The first time this remote reference is used, Singularity may display informational messages while it retrieves and prepares the image:

```text
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
[...]
INFO:    Creating SIF file...
```
{: .output}

Singularity then starts a container and executes the runscript provided by the image.

When Singularity processes the command, it:

1. reads the online registry reference
2. retrieves the Docker/OCI image manifest and filesystem layers from Docker Hub
3. converts the image into the Singularity Image Format
4. stores the converted image in its internal cache
5. starts a container from the cached image
6. executes the runscript defined by the image publisher

Run the same command again:

```bash
$ singularity run docker://docker.io/sylabsio/lolcow:latest
```
{: .source}

The second execution should start sooner because Singularity can reuse the converted image stored in its internal cache. In both executions, the container runs locally. The online registry reference tells Singularity where to retrieve the image, but it does not mean that the image is executed remotely.

> ## Local SIF file or online image reference?
>
> The two commands use the same published image, but they manage it differently.
>
> Run the SIF file stored in your personal image library:
>
> ```bash
> $ singularity run "$SINGULARITY_IMAGE"
> ```
> {: .source}
>
> Run the image using its Docker Hub reference:
>
> ```bash
> $ singularity run docker://docker.io/sylabsio/lolcow:latest
> ```
> {: .source}
>
> Using an online registry reference is convenient for quickly testing an image. Singularity manages the converted image in its internal cache, where cached objects may be identified by content-based hashes rather than recognisable image names.
>
> For images that you intend to retain and use in research workflows, prefer an explicitly named SIF file stored in your personal or project image library.
{: .callout}

</details>

> ## Docker/OCI image, SIF image, Docker Hub, Docker, and Singularity
>
> These related concepts should not be confused:
>
> - A **Docker/OCI image** is an image packaged according to the Docker/OCI image format. In a registry, it is commonly distributed as a manifest, configuration metadata, and a set of filesystem layers.
> - A **SIF image** is a container image stored in the Singularity Image Format. It is normally represented by a single `.sif` file.
> - **Docker Hub**, identified by `docker.io`, is an online registry that stores and distributes Docker/OCI images.
> - **Docker** is a container platform and engine that can build, distribute, and run Docker/OCI images (not used in this episode).
> - **Singularity** is the container engine used in this episode to run images. It can retrieve a Docker/OCI image from Docker Hub or another compatible registry, convert it into SIF, and run the resulting container without using the Docker engine.
>
> In this example:
>
> ```text
> Docker Hub
>     |
>     |  Docker/OCI image:
>     |  docker.io/sylabsio/lolcow:latest
>     v
> Singularity retrieves and converts the image
>     |
>     |  Singularity SIF image:
>     |  lolcow_latest.sif
>     v
> Singularity runs the container
> ```
> {: .output}
>
> Docker Hub is the source registry, the Docker/OCI image is the source image, and `lolcow_latest.sif` is the converted image to a `.sif` file managed locally. Singularity performs the conversion and runs the container. The Docker engine is not involved.
{: .callout}

### Running a command in a container

`singularity run` executes the default action defined by the image publisher. To execute a command of your choice, use `singularity exec`:

```bash
$ singularity exec "$SINGULARITY_IMAGE" cowsay "Hello from my local SIF image!"
```
{: .source}

The general form is:

```bash
$ singularity exec IMAGE COMMAND [ARGUMENTS...]
```
{: .source}

The image is followed by the command to execute and any arguments to pass to it.

Ask the `cowsay` application for help:

```bash
$ singularity exec "$SINGULARITY_IMAGE" cowsay -h
```
{: .source}

List the figures packaged with `cowsay`:

```bash
$ singularity exec "$SINGULARITY_IMAGE" cowsay -l
```
{: .source}

> ## Choose a figure
>
> Select one of the figures reported by `cowsay -l` and print your own message. For example, if the image includes the `dragon` figure:
>
> ```bash
> $ singularity exec "$SINGULARITY_IMAGE" cowsay -f dragon "Running from a container!"
> ```
> {: .source}
{: .challenge}

> ## Running `cowsay` without a message
>
> If you run `cowsay` without providing a message:
>
> ```bash
> $ singularity exec "$SINGULARITY_IMAGE" cowsay
> ```
> {: .source}
>
> the program waits for input from the terminal. Type your message, press <kbd>Enter</kbd>, and then press <kbd>Ctrl</kbd>+<kbd>D</kbd> to indicate the end of the input. `cowsay` will then display the message.
>
> Pressing <kbd>Ctrl</kbd>+<kbd>C</kbd> interrupts and cancels the command instead.
>
> Providing the message as a command-line argument is usually simpler.
{: .callout}

### Inspecting the container environment

The image packages both its applications and the user-space environment required by them. Compare the operating-system information visible on the host with that inside the container.

On the host:

```bash
$ cat /etc/os-release
```
{: .source}

Inside the container:

```bash
$ singularity exec "$SINGULARITY_IMAGE" cat /etc/os-release
```
{: .source}

The outputs may describe different Linux distributions or releases. The information printed inside the container describes the user-space environment packaged in the image. The container does not boot its own kernel. Its processes continue to use the host system's Linux kernel, as discussed in the introductory episode.

You can also use `which` to locate the three commands used by the image's default action:

```bash
$ singularity exec "$SINGULARITY_IMAGE" which date cowsay lolcat
```
{: .source}

```text
/bin/date
/usr/games/cowsay
/usr/games/lolcat
```
{: .output}

These paths belong to the container's filesystem. The commands do not need to be installed on the host. This output is different from that obtained when trying to locate the programs in the host.

### Running shell expressions with `bash -c`

The command passed directly to `singularity exec` must be an executable that Singularity can start. Some useful shell operations are not separate executable files.

For example, if we invoque `command -v` directly as another mean to locate the important tools in the image, we would get an error:

```bash
$ singularity exec "$SINGULARITY_IMAGE" command -v date cowsay lolcat
```
{: .source}

This fails because `command` is a shell built-in that reports how a shell would resolve one or more command names. But it is not a separate executable that Singularity can start.

To use a shell built-in like this one, start Bash inside the container and use its `-c` option:

```bash
$ singularity exec "$SINGULARITY_IMAGE" bash -c 'command -v date cowsay lolcat'
```
{: .source}

```text
/bin/date
/usr/games/cowsay
/usr/games/lolcat
```
{: .output}

The `-c` option tells Bash to interpret the following quoted string as a shell command. The quotation marks keep the complete command string together so that it can be passed to Bash inside the container.

This pattern is useful whenever the operation to execute inside a container includes shell features such as:

- shell built-ins
- pipelines
- redirections
- internal variable expansion
- multiple commands

> ## Add colour to the cow
>
> The output produced by `cowsay` is not colourful, while the image's default action produces colourful output.
>
> The image also contains `lolcat`, which adds terminal colours to text. How can you pass the output from `cowsay` to `lolcat` to end with a colourful message?
>
> > ## Naive solution
> >
> > A first attempt may be:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" cowsay "Hello from my local SIF image!" | lolcat
> > ```
> > {: .source}
> >
> > This does not run the complete pipeline inside the container. The host shell interprets the pipe before Singularity starts:
> >
> > ```text
> > singularity exec "$SINGULARITY_IMAGE" cowsay "Hello..."  |  lolcat
> >             runs inside the container                       runs on the host
> > ```
> > {: .output}
> >
> > Therefore, `cowsay` runs inside the container, but the host shell searches for `lolcat` on the host. The command fails if `lolcat` is not installed in the host.
> {: .solution}
>
> > ## Solution
> >
> > Pass the complete pipeline as a quoted command string to Bash inside the container:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" bash -c 'cowsay "Hello from my local SIF image!" | lolcat'
> > ```
> > {: .source}
> >
> > The host shell passes the quoted pipeline as one argument to `bash -c`. Bash inside the container interprets the pipe, so both `cowsay` and `lolcat` run inside the container.
> {: .solution}
{: .challenge}

> ## Reproduce the image's default action
>
> The image's default action uses:
>
> - `date` to print the current date and time
> - `cowsay` to place that text in an ASCII-art speech bubble
> - `lolcat` to add terminal colours
>
> Use `singularity exec` and `bash -c` to combine these commands into a pipeline that reproduces the image's default action.
>
> > ## Solution
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" bash -c 'date | cowsay | lolcat'
> >
> > ```
> > {: .source}
> >
> > Bash runs inside the container and interprets both pipe operators. Therefore, all three commands are resolved and executed using the container environment.
> >
> > This pipeline produces the same type of output as:
> >
> > ```bash
> > $ singularity run "$SINGULARITY_IMAGE"
> > ```
> > {: .source}
> >
> > `singularity run` executes the pipeline already defined in the image's runscript, while `singularity exec` with `bash -c` specifies the pipeline explicitly.
> {: .solution}
{: .challenge}

### Opening an interactive shell in a container

The `singularity exec` command runs a specified command and then returns control to the host shell. For interactive inspection and troubleshooting, use `singularity shell` instead:

```bash
$ singularity shell "$SINGULARITY_IMAGE"
```
{: .source}

The prompt changes to indicate that you are interacting with a shell inside the container environment:

```text
Singularity>
```
{: .output}

Commands entered at this prompt are interpreted by the shell running inside the container. For example, inspect the packaged operating-system environment:

```bash
Singularity> cat /etc/os-release
```
{: .source}

Because a shell is already running inside the container, shell built-ins such as `command` can be used directly without `bash -c`:

```bash
Singularity> command -v date cowsay lolcat
```
{: .source}

```text
/bin/date
/usr/games/cowsay
/usr/games/lolcat
```
{: .output}

Shell operators are also interpreted inside the container. Therefore, the pipeline used in the previous challenge can be entered directly:

```bash
Singularity> date | cowsay | lolcat
```
{: .source}

Similarly, you can provide your own colourful message:

```bash
Singularity> cowsay "Running interactively" | lolcat
```
{: .source}

The important difference is that the interactive container shell now interprets the commands and pipe operators. There is no need to start another shell with `bash -c`.

Exit the container shell when finished:

```bash
Singularity> exit
```
{: .source}

You can also press <kbd>Ctrl</kbd>+<kbd>D</kbd> to exit.

An interactive shell is useful for:

- inspecting the contents of an image
- locating applications and files
- testing commands
- investigating unexpected behaviour

For repeatable workflows and batch jobs, prefer `singularity exec`. Commands passed to `exec` can be recorded directly in scripts, while commands entered interactively are not automatically preserved.

A normal SIF image is read-only during execution. Exploring the image or attempting to modify its packaged files from an interactive shell does not permanently change the original SIF image.

### Image tags and reproducibility

An image tag identifies a published image variant. The `latest` tag is only a conventional name. It does not guarantee that an image contains the newest application version, and its contents may change when the publisher updates it.

The `sylabsio/lolcow` image is used with `latest` because that is the tag provided for this teaching example. For research and production workflows, prefer a meaningful version tag when one is available:

```text
docker://docker.io/namespace/application:1.2.3
```
{: .output}

A versioned tag communicates the intended software version more clearly, although publishers can technically update tags. Retaining the downloaded SIF file preserves the exact image contents that you obtained at that time.

### Getting help

Use `singularity help` to display general help:

```bash
$ singularity help
```
{: .source}

Add a command name for command-specific help:

```bash
$ singularity help pull
$ singularity help run
$ singularity help exec
$ singularity help shell
```
{: .source}

These commands display help for the Singularity container engine. To display help for an application packaged inside an image, execute that application's help command through the container, for example:

```bash
$ singularity exec "$SINGULARITY_IMAGE" cowsay -h
```
{: .source}

> ## Review the basic workflow
>
> In this episode, you followed the basic lifecycle of an existing container image:
>
> 1. Find an image in a registry.
> 2. Inspect its publisher, tags, documentation, and compatibility.
> 3. Pull it into a named SIF file for regular use.
> 4. Use `run`, `exec`, and `shell` for different interactions with the image.
> 5. Use `bash -c` when an operation requires shell syntax.
> 6. Optionally, use an online registry reference for a quick test.
>
> For your own workflows, store reusable SIF files in an organised location and record the original registry reference and tag from which each file was obtained.
{: .callout}
