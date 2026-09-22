---
title: "Building Docker/OCI images for HPC with Singularity"
teaching: 30
exercises: 20
questions:
- "Why do we build images with Docker but run them with Singularity on a shared HPC system?"
- "How can I describe and build a reproducible container image with a Dockerfile?"
- "How can I test the image locally before using it on a cluster?"
- "How can I publish an image to a registry and convert it to a SIF file?"
objectives:
- "Explain why Docker and Singularity have different roles in an HPC container workflow"
- "Read and write a basic Dockerfile using `FROM`, `LABEL`, `RUN`, `ENV`, `COPY`, and `CMD`"
- "Build and test an `amd64` Docker/OCI image"
- "Apply basic practices for build contexts, package installation, image tags, and runtime users"
- "Build an MPI application image from a Pawsey-provided base image"
- "Publish an image to Docker Hub and pull it as a named SIF file on an HPC system"
keypoints:
- "Docker is commonly used to build and test Docker/OCI images on a developer-controlled computer, while Singularity provides the unprivileged runtime and host integration required on a shared HPC system"
- "A Dockerfile records the base image and the instructions used to assemble a new image"
- "`COPY` adds files from the build context to the image"
- "Docker image layers support build caching, but Dockerfile instruction order and cleanup affect build efficiency and image size"
- "Use a small build context, a `.dockerignore` file, a trusted base image, and a meaningful image tag"
- "Build for the CPU architecture of the target HPC system"
- "Pawsey-provided base images offer tested starting environments for MPI and GPU applications on Pawsey systems"
- "Do not store passwords, access tokens, or other secrets in a Dockerfile, build argument, or image layer"
- "A registry provides the normal bridge between a Docker build environment and Singularity on an HPC system"
- "Singularity can pull a Docker/OCI image from a registry and convert it into a named SIF file"
---

### Prepare for the hands-on exercise

> ## Run on your local computer
>
> Run the following commands in the terminal on the local computer where Docker was installed and tested. Do not run these Docker commands on a Setonix login or compute node. This location remains in effect until another location callout appears.
{: .callout}

In your own computer, go to an adecuate directory and clone the training repository if you have not already done so:

```bash
$ mkdir -p "$HOME/containersTrainingPawsey"
$ cd "$HOME/containersTrainingPawsey"
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
```
{: .source}

Now move to the working directory for this episode:

```bash
$ cd "${TUTO}/demos/build_lolcow_docker"
$ pwd
```
{: .source}

The working directory should end with:

```text
singularity-containers/demos/build_lolcow_docker
```
{: .output}

### Use of Docker and Singularity

In the introductory episodes, we used Singularity to find, pull, and run existing Docker/OCI images on an HPC system. We also saw that Singularity converts registry images into read-only SIF files and runs container processes with the user's normal identity.

This episode introduces the other side of the workflow: creating a Docker/OCI image with Docker on a workstation, testing it, publishing it to a registry, and then pulling it with Singularity on the cluster.

The workflow is:

```text
Dockerfile
    |
    |  docker build
    v
Local Docker/OCI image
    |
    |  docker push
    v
Container registry
    |
    |  singularity pull
    v
Named SIF image on the HPC system
```
{: .output}

#### Why build with Docker for an HPC workflow?

Singularity is the container engine used to run images throughout this training because it is designed for shared HPC systems and integrates containerised applications with host filesystems, schedulers, MPI libraries, GPUs, and high-speed interconnects.

For building images, this training uses Docker and Dockerfiles. The first important advantage is Docker's layered build model. Docker records filesystem changes from build instructions in image layers and can reuse unchanged layers from its build cache. During development, editing a later Dockerfile instruction may therefore require rebuilding only that instruction and the instructions that follow it, rather than repeating the complete build. This makes the iterative cycle of editing, building, and testing more efficient.

Singularity also caches downloaded Docker/OCI layers and converted images, so repeated pulls do not necessarily download the same content again. However, that is different from Docker's instruction-level build cache: a Singularity definition-file build does not provide the same Dockerfile layer-by-layer workflow for incrementally rebuilding a customised image.

The second important advantage is interoperability. Docker builds images in the widely supported Docker/OCI ecosystem. These images can be stored in standard OCI-compatible registries and used by many container tools. Depending on the tool, an image may be run directly, imported, or converted into its native format. Singularity, for example, can retrieve a Docker/OCI image from a registry, assemble its layers, and convert it into a SIF image for execution on an HPC system.

In this sense, Docker/OCI images are a broadly interoperable distribution format. This does not mean that every image behaves identically with every container engine. Runtime features, image metadata, security models, CPU architecture, and host integration can differ. The final image must still be tested with Singularity on the target HPC system.

The different security models also matter. A traditional Docker installation uses a daemon that normally runs with root privileges and can perform host-level operations such as starting containers and mounting host directories. Giving general users control of that daemon would provide highly privileged access on a shared machine. Docker is therefore suitable for a participant's own computer, but it is not provided as the general user-facing container runtime on Setonix.

Singularity is designed for shared HPC systems. Normal container execution uses the invoking user's host identity rather than requiring the user to control a privileged daemon, and the Pawsey modules integrate containers with the filesystems, scheduler, MPI libraries, GPUs and other host facilities. This is why the training uses two engines rather than running Docker directly on Setonix.

The tools therefore have complementary roles in this training:
- **Docker builds and tests the Docker/OCI image** on the participant's local computer, where the participant controls the development environment.
- **A container registry stores and distributes the image** in a widely supported format.
- **Singularity retrieves and converts the image to SIF, then runs it** using the shared-system execution model on Setonix.

The registry transfer and conversion to SIF add steps, but they preserve Docker's mature Dockerfile and build-cache workflow without requiring Docker as the cluster runtime. They also create an explicit, immutable SIF artefact that can be organised and validated on the target system.

Docker is not used to run the workload on Setonix. A successful local Docker test cannot establish that the image works with Singularity, the user's Setonix identity, read-only SIF execution, MPI, GPUs or Pawsey's host integration. The final image must therefore be validated with Singularity and the supported launch model on Setonix.

### Prepare the Docker build directory

Perform the Docker sections of this episode on the computer where Docker was installed and tested in the installation episode. Do not run these commands on a Setonix login or compute node.

List the files in the directory:

```bash
$ ls -l
```
{: .source}

The relevant entries should look like:

```text
Dockerfile -> lolcow.dockerfile
lolcow.dockerfile
lolcow-message.txt
```
{: .output}

The actual recipe is named `lolcow.dockerfile` so that its purpose remains clear when it is viewed outside this directory or alongside recipes for other images. The symbolic link named `Dockerfile` points to that recipe. `Dockerfile` is the default filename used by Docker, so the link allows the standard `docker build ... .` command to find the recipe without an additional option.

The `lolcow-message.txt` file contains the message that will be copied into the image and displayed by the default container action. Near the end of the build, the recipe and message file are also preserved inside the image for later inspection.

### Read the Dockerfile

The `lolcow.dockerfile` recipe is:

```dockerfile
# Start from a versioned Ubuntu image on Docker Hub
FROM docker.io/ubuntu:24.04

# Install the applications and remove package-manager cache files
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Make the installed commands available by name at runtime
ENV PATH="/usr/games:${PATH}"

# Copy the message displayed by the default container action
COPY lolcow-message.txt /usr/local/share/lolcow/message.txt

# Preserve the recipe and build input files inside the image
ARG IMAGE_BUILD_INFO_DIR="/opt/build-info-and-recipes/lolcow"
RUN mkdir -p "${IMAGE_BUILD_INFO_DIR}"
COPY lolcow.dockerfile \
     lolcow-message.txt \
     "${IMAGE_BUILD_INFO_DIR}/"

# Record standard OCI image metadata and the build-information location
LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre" \
      org.opencontainers.image.source="https://github.com/PawseySC/singularity-containers" \
      au.org.pawsey.image.build-info-dir="${IMAGE_BUILD_INFO_DIR}"


# Define the default action for docker run and singularity run
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]
```
{: .source}

Docker reads the instructions from top to bottom. Each instruction describes part of the resulting image or its default runtime configuration.

#### `FROM`: select a base image

```dockerfile
FROM docker.io/ubuntu:24.04
```
{: .source}

`FROM` begins a build stage and selects its base image. This example starts from the versioned `docker.io/ubuntu:24.04` image rather than `docker.io/ubuntu:latest`. The explicit `docker.io` component identifies Docker Hub as the registry.

Choose base images from trusted publishers and prefer a supported, suitably small image that provides what the application needs. A versioned tag communicates the intended base more clearly, although tags can still be updated by their publisher. For stricter provenance, production builds may pin the base image by digest and update that digest deliberately.

#### `RUN`: execute build-time commands

```dockerfile
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*
```
{: .source}

`RUN` executes commands while the image is being built and stores the resulting filesystem changes in an image layer.

This instruction:

1. enables strict and verbose shell behaviour for the build step
2. updates the Ubuntu package index
3. installs the two required packages without additional recommended packages
4. removes cached package data that is not needed at runtime

The package-index update, installation, and cleanup are performed in the same `RUN` instruction. Removing files in a later layer would not remove them from the earlier layer in which they were created.

Using fewer `RUN` instructions does not by itself guarantee a good Dockerfile. Combine commands when their changes belong in one filesystem layer, but keep the result readable and ensure failures stop the build.

##### Combining commands in one `RUN` instruction

The package-index update, installation, and cleanup are kept in one `RUN` instruction so their final filesystem changes are recorded in one image layer. This prevents temporary package data from remaining in an earlier layer.

In this recipe:

- `\` continues the Dockerfile instruction onto the next physical line; it does not join shell commands by itself
- `;` separates one shell command from the next
- `set -e` stops the shell after an unhandled command failure, which is important because semicolon-separated commands would otherwise continue
- `set -u` reports an error when an unset variable is expanded
- `set -x` prints expanded commands as they execute, making the build log easier to diagnose

Together, `set -eux` provides failure handling, unset-variable checking, and command tracing for the remainder of this `RUN` instruction. Because `-x` prints expanded commands, do not use traced commands to handle passwords, access tokens, or other secrets.

> ## Optional: Using `&&` instead
>
> A more common Dockerfile style connects the same commands with `&&`:
>
> ```dockerfile
> RUN export DEBIAN_FRONTEND=noninteractive \
>     && apt-get update \
>     && apt-get install -y --no-install-recommends \
>         cowsay \
>         lolcat \
>     && apt-get clean \
>     && rm -rf /var/lib/apt/lists/*
> ```
> {: .source}
>
> The `&&` operator runs the next command only when the preceding command succeeds. For this simple linear sequence, it provides explicit failure chaining without requiring `set -e`.
>
> Both styles create one image layer because each uses one `RUN` instruction. This training keeps the `set -eux` style because it is also suitable for longer instructions that use variables, loops, or conditionals.
{: .solution}

#### `ENV`: define a runtime environment variable

```dockerfile
ENV PATH="/usr/games:${PATH}"
```
{: .source}

`ENV` defines an environment variable that persists in the image and is normally present when a container is started. Ubuntu installs `cowsay` under `/usr/games`, so this instruction adds that directory to `PATH`.

An `export` performed inside one `RUN` instruction affects only the shell used for that build step. Use `ENV` when the setting should form part of the image's runtime environment.


#### `COPY`: add a file from the build context

```dockerfile
COPY lolcow-message.txt /usr/local/share/lolcow/message.txt
```
{: .source}

`COPY` adds files or directories from the build context to the image. Here, the source is `lolcow-message.txt` in the Docker build directory, and the destination is `/usr/local/share/lolcow/message.txt` inside the image. Docker creates the required destination directories when it performs the copy.

The source must be present in the build context and must not be excluded by `.dockerignore`. Unlike a bind mount, the copied file becomes part of the built image and remains available when the image is transferred to another system.

Another instruction that could copy this local file is `ADD`. However, `ADD` has additional capabilities, such as automatically extracting local tar archives and retrieving remote sources. For straightforward copies of local files and directories, prefer `COPY` because its behaviour and intent are clearer.

#### Preserve the recipe, build inputs, and metadata

Near the end of the recipe, `ARG`, `RUN`, and `COPY` preserve the Dockerfile and the files consumed by the build under `/opt/build-info-and-recipes/lolcow`. The final `LABEL` instruction records standard OCI metadata, the source repository, and the internal build-information directory.

Keeping these frequently edited instructions near the end allows Docker to reuse the earlier package-installation layers when the preserved recipe or labels change. The embedded files support later inspection but do not replace the version-controlled source repository, image digest, base-image digest, and recorded build command.

#### `CMD`: define the default action

```dockerfile
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]
```
{: .source}

`CMD` is the final instruction in this recipe and defines the default command used when a Docker container is started without another command. Placing it last makes the final runtime decision easy to locate. The JSON, or exec, form preserves the command and arguments as separate values.

The redirection and pipeline must be interpreted by a shell, so the Dockerfile explicitly starts `bash -c`. `cowsay` reads the copied message through standard input, and `lolcat` processes the resulting output. This is the same shell-expression principle used with `singularity exec` in the basic Singularity episode.

Later in the episode, the `CMD` instruction is updated to use `lolcat --force`, ensuring that colour codes are emitted even when Docker has not allocated a terminal.

> ## `CMD` is a default, not a build command
>
> `RUN` executes while the image is built. `CMD` records the default command to execute later when a container is started. A Dockerfile can contain only one effective `CMD`; if several are present, the last one takes effect.
{: .callout}

### Build the image for the target HPC architecture

The general form of the build command is:

```text
docker build [OPTIONS] BUILD_CONTEXT
```
{: .output}

The command builds a Docker/OCI image by processing a Dockerfile and the files available in the build context:

- `[OPTIONS]` configures the build, including the target platform, image reference, and Dockerfile to use.
- `BUILD_CONTEXT` identifies the directory, URL, or standard-input stream containing the files available to the build. A local directory is the most common choice.

The build context is a required positional argument and normally appears last. Dockerfile instructions such as `COPY` and `ADD` can access files from this context. Docker normally looks for a file named `Dockerfile` at the root of the context unless another recipe is selected with `--file`.

Define the complete local image reference, including its repository name and tag:

```bash
$ COW_IMAGE="lolcow:2026.09"
```
{: .source}

Setonix compute nodes use the `amd64` architecture, also called `x86_64`. Build explicitly for that target:

```bash
$ docker build --platform linux/amd64 --tag "$COW_IMAGE" .
```
{: .source}

In this command:

- `--platform linux/amd64` selects the target operating system and CPU architecture.
- `--tag "$COW_IMAGE"` assigns the repository name and tag `lolcow:2026.09` to the resulting local image. The short form of `--tag` is `-t`.
- `.` selects the current directory as the build context.

The `COW_IMAGE` variable contains the local repository name `lolcow` and tag `2026.09`. Docker finds the recipe through the symbolic link named `Dockerfile`, while the final `.` makes the files in the current directory available to instructions such as `COPY` and `ADD`. The context should therefore contain only what the build requires.

> ## Limiting a larger build context
>
> This example has a small, controlled build context containing only the recipe and its message file, so it does not require a `.dockerignore` file. Larger projects commonly use `.dockerignore` to reduce unnecessary build-context processing and to prevent broad `COPY` or `ADD` instructions from including unwanted or sensitive files.
{: .callout}

List the newly built image:

```bash
$ docker image ls "$COW_IMAGE"
```
{: .source}

The output should include the image name and tag:

```text
IMAGE            ID             DISK USAGE   CONTENT SIZE   EXTRA
lolcow:2026.09   fee6b76c45f3   53.4MB       53.4MB
```
{: .output}

The image ID and reported sizes may differ.

Confirm the operating system and CPU architecture recorded for the image:

<!-- The raw block prevents Jekyll/Liquid from interpreting Docker's Go-template braces. -->
```bash
$ docker image inspect "$COW_IMAGE" --format '{% raw %}{{.Os}}/{{.Architecture}}{% endraw %}'
```
{: .source}

The expected output is:

```text
linux/amd64
```
{: .output}

Now run the image's default action without `--rm`:

```bash
$ docker run "$COW_IMAGE"
```
{: .source}

The output should contain the message copied from `lolcow-message.txt`:

```text
 _________________________________________
/ Built with Docker and ready to run with \
\ Singularity!                            /
 -----------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
{: .output}

The exact spacing may vary. At this stage, the output may have no colours because `docker run` does not allocate a terminal by default and `lolcat` may suppress colour when its output is not connected to a terminal.

The `cowsay` pipeline has finished, so no process from this container is still running. Confirm that the default container list is empty:

```bash
$ docker container ls
```
{: .source}

Now include stopped containers in the listing:

```bash
$ docker container ls --all
```
{: .source}

The output should contain a stopped container created from the image. Docker assigns a generated name when `--name` is not specified:

```text
CONTAINER ID   IMAGE            COMMAND                  CREATED          STATUS                      PORTS   NAMES
a1b2c3d4e5f6   lolcow:2026.09   "bash -c 'cowsay …'"   10 seconds ago   Exited (0) 8 seconds ago           generated_name
```
{: .output}

The container ID, generated name, and times will differ.

`docker run` creates a new container from the image and starts its main process. When that process finishes, the container stops, but Docker retains the container object by default. The retained object includes its configuration, logs, metadata, and writable filesystem layer. This allows a stopped container to be inspected, restarted, or used to recover files created during its execution.

For these short tests, there is no useful state to retain. Remove all exited containers:

```bash
$ docker container rm $(docker container ls --all --quiet --filter status=exited)
```
{: .source}

The command substitution is intentionally unquoted so that each container ID is passed to `docker container rm` as a separate argument.

Run the image again, this time with automatic cleanup:

```bash
$ docker run --rm "$COW_IMAGE"
```
{: .source}

The `--rm` option instructs Docker to remove the container automatically after its main process exits. Confirm that this second test did not leave another stopped container:

```bash
$ docker container ls --all
```
{: .source}

For short-lived tests in this episode, continue using `--rm` unless you deliberately need to inspect or restart the stopped container afterward.

> ## Running the `amd64` image on an Apple Silicon Mac
>
> On an `arm64` computer, such as an Apple Silicon Mac, Docker may report that the requested image platform does not match the detected host platform. This is expected because the image was deliberately built for the `linux/amd64` architecture used on Setonix. Docker Desktop can use emulation to run it.
>
> Make the intended runtime platform explicit with:
>
> ```bash
> $ docker run --rm --platform linux/amd64 "$COW_IMAGE"
> ```
> {: .source}
>
> This suppresses the platform-mismatch warning but still uses emulation on an `arm64` host.
{: .callout}


### Alternative build commands

The explicit `linux/amd64` build used above is the recommended command for preparing this image for Setonix. The following alternatives illustrate what happens when the target platform is omitted, how Docker selects its default recipe, and how Docker Buildx can validate a build configuration.

#### Build for the local computer's native architecture

If `--platform` is omitted, Docker normally builds for the builder's native platform. Use a different image tag so that this comparison does not replace the `linux/amd64` image required for Setonix:

```bash
$ NATIVE_IMAGE="lolcow:native"
$ docker build --tag "$NATIVE_IMAGE" .
```
{: .source}

Inspect the resulting image:

<!-- The raw block prevents Jekyll/Liquid from interpreting Docker's Go-template braces. -->
```bash
$ docker image inspect "$NATIVE_IMAGE" --format '{% raw %}{{.Os}}/{{.Architecture}}{% endraw %}'
```
{: .source}

The result normally reflects the local builder:

- `linux/amd64` on a typical x86-64 Linux or Windows system
- `linux/arm64` on an Apple Silicon Mac

A native `arm64` image may be convenient for local execution, but it is not the image that this exercise prepares for Setonix.

#### Select the recipe explicitly

By default, `docker build` looks for a file named `Dockerfile` at the root of the build context. In this example, that conventional name is a symbolic link to the descriptively named recipe `lolcow.dockerfile`.

Remove the symbolic link:

```bash
$ rm Dockerfile
```
{: .source}

Try the default build command again:

```bash
$ docker build --platform linux/amd64 --tag "$COW_IMAGE" .
```
{: .source}

The build should fail with an error similar to:

```text
failed to read dockerfile: open Dockerfile: no such file or directory
```
{: .error}

The exact error may differ between Docker versions. The build fails because no file named `Dockerfile` is now present at the root of the build context.

Select the actual recipe explicitly with `--file`:

```bash
$ docker build --platform linux/amd64 \
    --file lolcow.dockerfile --tag "$COW_IMAGE" .
```
{: .source}

Recreate the relative symbolic link and inspect it:

```bash
$ ln -s lolcow.dockerfile Dockerfile
$ ls -l Dockerfile
```
{: .source}

The output should show:

```text
Dockerfile -> lolcow.dockerfile
```
{: .output}

Using a descriptive filename such as `lolcow.dockerfile` makes the recipe identifiable when it is viewed outside its original directory or alongside other recipes. Otherwise, a user can accumulate many unrelated files all named `Dockerfile`. The relative symbolic link preserves Docker's conventional default filename while keeping the actual recipe descriptive and portable with the repository.

#### Check the build configuration with Docker Buildx

Docker Buildx is a Docker CLI plugin that exposes extended capabilities of the BuildKit backend. It became widely associated with cross-platform and multi-platform builds because it provided direct access to configurable builders and the `--platform` option. In current Docker installations, the ordinary `docker build` command also uses BuildKit and can build directly for `linux/amd64`, as demonstrated earlier.

Buildx remains useful for additional operations such as build checks, multi-platform builds, configurable builders, cache import and export, and explicit output selection. For this episode, use its `--check` option to analyse the Dockerfile and build options without executing the complete image build:

```bash
$ docker buildx build --check --file lolcow.dockerfile .
```
{: .source}

Build checks act like Dockerfile linting: they report recognised issues, outdated practices, or inconsistencies in the build configuration. A clean recipe may complete without warnings. A reported warning does not necessarily mean that a normal build would fail.

The `--check` option requires Docker Buildx 0.15.0 or later. Check the installed version with:

```bash
$ docker buildx version
```
{: .source}

The main `docker build --platform linux/amd64 ...` command remains the recommended way to build this image during the episode. Buildx is introduced here so that you recognise the extended build interface and one of its practical development tools.

### Modify the image and reuse the build cache

Docker's layered build model allows unchanged build results to be reused. To demonstrate this, we will make a small change to the final `CMD` instruction, rebuild the image, and inspect which earlier steps Docker retrieves from its build cache. The practical change adds `--force` because `lolcat` normally suppresses colour when its standard output is not connected to a terminal.

Open `lolcow.dockerfile` in a text editor and change the final instruction from:

```dockerfile
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]
```
{: .source}

to:

```dockerfile
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat --force"]
```
{: .source}

The `--force` option tells `lolcat` to emit colour codes even when its output is not connected to a terminal.

Rebuild the image using the same name and tag:

```bash
$ docker build --platform linux/amd64 -t "$COW_IMAGE" .
```
{: .source}

Which build steps are reused from the cache, and which part changes?

> ## Solution
>
> The base image, package-installation layer, environment setting, and runtime copy of the message are unchanged, so Docker can reuse their cached results. Because the Dockerfile itself is preserved inside the image, changing its final `CMD` also invalidates the build-information `COPY`. The following `LABEL` and final `CMD` instructions are therefore processed again, while the earlier package-installation layers remain cached.
>
> This demonstrates why instructions that change frequently are normally placed after stable and expensive build steps.
{: .solution}

Run the rebuilt image:

```bash
$ docker run --rm "$COW_IMAGE"
```
{: .source}

The same message should now be displayed in colour. The terminal must support ANSI colour sequences, as standard terminals on current Linux, macOS, and Windows installations normally do.

### Run other commands with Docker

Override the default `CMD` by supplying another command after the image name:

```bash
$ docker run --rm "$COW_IMAGE" cat /etc/os-release
```
{: .source}

Docker runs the selected command instead of the image's default `CMD`.

### Open an interactive shell with Docker

As with `singularity shell`, Docker can open an interactive shell for inspecting and testing an image:

```bash
$ docker run --rm --interactive --tty "$COW_IMAGE" bash
```
{: .source}

The `--interactive` option keeps standard input open, while `--tty` allocates a pseudo-terminal. The final `bash` overrides this image's default `CMD` and starts an interactive Bash shell instead. Other images may already define a shell as their default action, may provide a different shell such as `sh`, or may require an `ENTRYPOINT` to be overridden.

The prompt should change to something similar to:

```text
root@CONTAINER-ID:/#
```
{: .output}

The container ID and exact prompt will differ.

Unlike normal Singularity execution on the cluster, this Docker container runs as root inside the container because the Dockerfile does not define another runtime user. This root identity applies within Docker's container environment. It is not the root user of the host operating system.

An image that works as Docker root can still fail on Setonix, where Singularity normally runs the process with your host user identity and the SIF filesystem is read-only. Install software into the image at build time, but design the application to write results, caches, temporary files, and runtime configuration only to writable locations such as the current working directory, `/tmp`, or bind-mounted project and scratch directories. The runtime must not require `sudo`, creation of system accounts, or modification of directories such as `/usr` and `/opt`. Detailed non-root and arbitrary-user testing is covered in the advanced Docker episode.

Inspect the packaged operating-system environment and locate the installed commands:

```bash
root@CONTAINER-ID:/# cat /etc/os-release
root@CONTAINER-ID:/# command -v cowsay lolcat
```
{: .source}

Run the installed applications directly:

```bash
root@CONTAINER-ID:/# cowsay "Running interactively with Docker"
root@CONTAINER-ID:/# cowsay "Running interactively with Docker" | lolcat --force
```
{: .source}

As in the Singularity episode, the interactive shell is useful for inspection and testing. Changes entered interactively are not recorded in the Dockerfile and should not replace a reproducible build.

### Access the host working directory

Singularity normally makes the host current working directory available inside the container automatically. Docker does not do this by default.

From the interactive Docker shell, inspect the current directory:

```bash
root@CONTAINER-ID:/# pwd
root@CONTAINER-ID:/# ls
```
{: .source}

The files from the host `build_lolcow_docker` directory are not visible. Exit the container:

```bash
root@CONTAINER-ID:/# exit
```
{: .source}

Start another interactive container and explicitly bind mount the host current directory at `/work`:

```bash
$ docker run \
    --rm \
    --interactive \
    --tty \
    --mount type=bind,source="$PWD",target=/work \
    "$COW_IMAGE" \
    bash
```
{: .source}

Inside the container, inspect the mounted directory:

```bash
root@CONTAINER-ID:/# ls -l /work
```
{: .source}

The directory should contain files including:

```text
Dockerfile
lolcow.dockerfile
lolcow-message.txt
```
{: .output}

The files remain stored on the host. Docker only makes the host directory accessible at `/work` for this container.

Exit when finished:

```bash
root@CONTAINER-ID:/# exit
```
{: .source}

#### Copy a file from the image to the host non-interactively

As in the basic Singularity episode, a command can copy a file packaged inside the image to a host directory without opening an interactive shell. Docker does not mount the host current working directory automatically, so make it available at `/work` and copy the packaged message into it:

```bash
$ docker run \
    --rm \
    --mount type=bind,source="$PWD",target=/work \
    "$COW_IMAGE" \
    cp /usr/local/share/lolcow/message.txt /work/lolcow-message.copy.txt
```
{: .source}

The source path is the file copied into the image during the build. The destination is within `/work`, which maps to the host current working directory.

After the container exits, inspect the copied file from the host:

```bash
$ cat lolcow-message.copy.txt
```
{: .source}

The output should be:

```text
Built with Docker and ready to run with Singularity!
```
{: .output}

The container was removed automatically because `--rm` was used, but the copied file remains because it was written through the bind mount to the host filesystem.

### Use Pawsey-provided base images

Pawsey publishes container base images that users can extend for their own applications. These are tested starting environments for Pawsey systems, including MPICH-based images prepared for the hybrid MPI model on Setonix and ROCm-based images for AMD GPU workloads. Each derived application image and workflow must still be validated on the target system.

The recipes are available from the [Pawsey container recipes repository](https://github.com/PawseySC/pawsey-containers), while the corresponding Docker/OCI images are published under the [Pawsey organisation on Quay.io](https://quay.io/pawsey).

Move to the separate Mandelbrot build context:

```bash
$ cd "$TUTO/demos/build_mandelbrot_docker"
$ ls -l
```
{: .source}

The directory contains:

```text
README.md
THIRD_PARTY_NOTICES.md
mandelbrot_mpi.dockerfile
mpi-mandelbrot.cpp
runMandelbrotSingularityPawsey.slurm.sh
runMandelbrotDocker.sh
```
{: .output}

The image is built from:

```dockerfile
FROM quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04
```
{: .source}

The base image supplies MPICH, compiler wrappers, and MPI runtime tools. The derived image compiles `mpi-mandelbrot.cpp` with `mpic++`, installs ImageMagick for PPM-to-PNG conversion, and preserves `THIRD_PARTY_NOTICES.md` inside the image. The recipe and all image-build inputs are also preserved under `/opt/build-info-and-recipes/mandelbrot-mpi`.

The complete recipe is:

```dockerfile
# Build the application on Pawsey's Setonix-compatible MPICH base image
FROM quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04

# Install the utility used to convert the PPM result to PNG
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends imagemagick; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Copy and compile the MPI application with the compiler from the base image
COPY mpi-mandelbrot.cpp /tmp/mpi-mandelbrot.cpp
RUN mpic++ \
        -std=c++17 \
        -O3 \
        -Wall \
        -Wextra \
        -Wpedantic \
        -o /usr/local/bin/mpi-mandelbrot \
        /tmp/mpi-mandelbrot.cpp \
    && rm -f /tmp/mpi-mandelbrot.cpp

# Preserve third-party acknowledgements and licence information
COPY THIRD_PARTY_NOTICES.md \
    /usr/local/share/doc/mpi-mandelbrot/THIRD_PARTY_NOTICES.md

# Preserve the recipe and build input files inside the image
ARG IMAGE_BUILD_INFO_DIR="/opt/build-info-and-recipes/mandelbrot-mpi"
RUN mkdir -p "${IMAGE_BUILD_INFO_DIR}"
COPY mandelbrot_mpi.dockerfile \
     mpi-mandelbrot.cpp \
     THIRD_PARTY_NOTICES.md \
     "${IMAGE_BUILD_INFO_DIR}/"

# Record standard OCI image metadata and the build-information location
LABEL org.opencontainers.image.title="MPI Mandelbrot renderer" \
      org.opencontainers.image.description="MPI training application built on Pawsey's MPICH base image" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/PawseySC/singularity-containers" \
      au.org.pawsey.image.build-info-dir="${IMAGE_BUILD_INFO_DIR}"


# Display application help when no other command is supplied
CMD ["mpi-mandelbrot", "--help"]
```
{: .source}

The Mandelbrot Dockerfile and all files consumed by `docker build` are preserved under `/opt/build-info-and-recipes/mandelbrot-mpi`; the final labels record this location. The host-side run scripts are not copied because they are runtime launch policy, not image-build inputs.

The image deliberately contains the application rather than a site-specific MPI launcher. A simple wrapper could have been copied into the image to run `mpiexec`, execute the renderer, and convert its output. That approach is convenient for a self-contained Docker demonstration, but it does not match the supported Setonix launch model, where Slurm starts the tasks outside the image. Embedding the launcher would also make it part of the immutable image, so refining the launch policy would require rebuilding and redistributing the image.

Instead, this example keeps launch scripts on the host. The application image remains reusable, while each environment selects the correct launcher:

```text
Local Docker: host script -> docker run -> container mpiexec -> MPI ranks
Setonix:      Slurm srun -> one singularity exec per task -> MPI ranks
```
{: .output}

This separation is intentional: the image contains the application and stable dependencies; command-line arguments and environment variables configure a run; bind mounts provide input data and writable output locations; and the host-side script records environment-specific launch policy. Keeping data, credentials, site paths, and Slurm resource requests outside the image makes the image easier to reuse and validate.

#### Build the MPI application image

```bash
$ MPI_IMAGE="mandelbrot-mpi:2026.09"
$ docker build \
    --platform linux/amd64 \
    --file mandelbrot_mpi.dockerfile \
    --tag "$MPI_IMAGE" \
    .
```
{: .source}

Check the default application help:

```bash
$ docker run --rm --platform linux/amd64 "$MPI_IMAGE"
```
{: .source}

#### Test MPI locally with Docker

The host-side `runMandelbrotDocker.sh` script uses `mpiexec` inside one Docker container. This is a local functional test on one computer, not the Setonix launch method. Its default workload is `1200 x 800` pixels, 500 maximum iterations, and four MPI processes. `WIDTH`, `HEIGHT`, `ITERATIONS`, `MPI_PROCESSES`, `CENTRE_REAL`, `CENTRE_IMAGINARY`, and `SCALE` can be overridden as environment variables.

```bash
$ ./runMandelbrotDocker.sh
```
{: .source}

For example, render a zoomed view with a larger workload without editing the script:

```bash
$ MPI_PROCESSES=8 \
    WIDTH=3000 \
    HEIGHT=2000 \
    ITERATIONS=1000 \
    CENTRE_REAL=-0.743643887037151 \
    CENTRE_IMAGINARY=0.131825904205330 \
    SCALE=0.002 \
    ./runMandelbrotDocker.sh
```
{: .source}

The script defines `$OUTPUT_DIR` on the host and bind mounts it at `/output` inside Docker. Commands inside the container use `/output/$FILE_PPM` and `/output/$FILE_PNG`, while cleanup and reporting on the host use paths under `$OUTPUT_DIR`. The script starts four MPI ranks inside one container, converts the PPM result to PNG with a second serial container command, and removes the intermediate PPM file. Representative output is:

```text
MPI Mandelbrot renderer
Image size: 3000 x 2000
Maximum iterations: 1000
Centre: (-0.743644, 0.131826)
Scale: 0.002
MPI processes: 8
PPM output: /output/mandelbrot.docker.ppm
Rendering completed in 0.420 seconds
Created /path/to/build_mandelbrot_docker/output/mandelbrot.docker.png
```
{: .output}

Confirm and open the result:

```bash
$ ls -lh output/mandelbrot.docker.png
```
{: .source}

### Review image-building practices

The lolcow and Mandelbrot examples demonstrate practices that should be retained in research and HPC container workflows:

- Start from an image maintained by a trusted project, vendor, or organisation, and prefer explicit application and base-image versions over `latest` where practical. In this episode, the examples use `docker.io/ubuntu:24.04` and `quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04`.
- Keep separate, focused build contexts for unrelated images. Keep each context small and use `.dockerignore` when needed to exclude unnecessary or sensitive files.
- Record the Dockerfile and related build files in version control. Preserve the recipe and image-build inputs inside the image when this supports provenance and later inspection. Record standard OCI metadata and preserve licence notices with redistributed software.
- Use `COPY` for ordinary file and directory copies. Use `ADD` only when its additional behaviour is specifically required.
- Install only required packages and remove package-manager caches in the same `RUN` instruction that installs them.
- Place stable and expensive dependency steps before frequently changing application files to improve build-cache reuse.
- Do not place passwords, private keys, access tokens, licence files, or other secrets in the Dockerfile, build arguments, environment variables, or copied build context. Use the build system's supported secret mechanism when a build must access protected resources.
- Build explicitly for the architecture of the target system.
- Keep site-specific launch policy outside an immutable application image when the target environments require different launch mechanisms.
- Rebuild and test images regularly so that base-image and package security updates are incorporated.
- Test locally, then validate the final SIF image with the supported runtime and launch model on the target HPC system, including MPI, GPU, filesystem, correctness, and performance behaviour where relevant.

For compiled applications, consider a **Docker multi-stage build**. A Dockerfile can use multiple `FROM` instructions, commonly naming a build stage with `AS`. A later runtime stage can use `COPY --from=<stage>` to copy only the compiled application and other required artefacts, leaving compilers, source files, and development packages out of the final image. The final stage must still provide every runtime library required by the copied application.

Image size is also reduced by choosing an appropriate base image, installing only required packages, using `--no-install-recommends`, removing package-manager caches in the same `RUN` instruction, keeping the build context small, and excluding unnecessary files with `.dockerignore`. Combining related commands can prevent temporary files from remaining in an earlier layer, but merely reducing the number of layers does not by itself guarantee a smaller or better image. A complete multi-stage and image-size exercise is covered in the advanced Docker episode.

### Publish the Mandelbrot image to Docker Hub

A registry is the normal way to move a Docker/OCI image from the build computer to the HPC system. This section uses Docker Hub, which was introduced in the setup episode.

The [Setup Docker on your computer]({% link _episodes/00-setup-docker.md %}) episode covered creating and verifying a Docker Hub account and testing image publishing. Use the same Docker ID here. If `DOCKER_ID` is not defined in the current shell, assign it now. Replace `<docker-id>` with your Docker ID and do not include the angle brackets:

```bash
$ DOCKER_ID="<docker-id>"
$ MPI_REMOTE_IMAGE="docker.io/${DOCKER_ID}/mandelbrot-mpi:2026.09"
```
{: .source}

> ## Windows PowerShell syntax
>
> In Windows PowerShell, assign the variables with:
>
> ```powershell
> PS> $DOCKER_ID = "<docker-id>"
> PS> $MPI_REMOTE_IMAGE = "docker.io/${DOCKER_ID}/mandelbrot-mpi:2026.09"
> ```
> {: .source}
>
> The `$` characters in `$DOCKER_ID` and `$MPI_REMOTE_IMAGE` are part of the PowerShell variable names and must be typed. The later Docker commands use the same quoted variable references in PowerShell, Bash, and Zsh.
{: .solution}

Before pushing, sign in to Docker Hub in a web browser and create a **public** repository named `mandelbrot-mpi` under your Docker ID, following the same repository-creation process used for `first-image` in the setup episode. The complete repository name will be `docker.io/<docker-id>/mandelbrot-mpi`. Do not include the angle brackets when substituting your Docker ID.

Authenticate from the Docker client:

```bash
$ docker login docker.io
```
{: .source}

Follow the authentication instructions shown by Docker. Do not put the password or access token directly in the Dockerfile or shell history.

Add the registry-qualified tag to the existing local image:

```bash
$ docker tag "$MPI_IMAGE" "$MPI_REMOTE_IMAGE"
```
{: .source}

The new tag is another name for the same local image. Confirm both names:

```bash
$ docker image ls
```
{: .source}

Push the registry-qualified image:

```bash
$ docker push "$MPI_REMOTE_IMAGE"
```
{: .source}

After the push completes, inspect the repository and tag in Docker Hub. For this training exercise, the repository must be readable from the cluster. Do not publish proprietary, confidential, export-controlled, licensed, or otherwise restricted software in a public repository.

> ## Image tags and immutable identity
>
> A tag such as `2026.09` is a readable reference, but a registry publisher can move a tag to different image content. A pushed image is also identified by a content digest beginning with `sha256:`.
>
> For a reproducible workflow, record the registry, namespace, repository, tag, digest, Dockerfile revision, and relevant build inputs. A retained SIF file also preserves the exact content that was pulled at that time.
{: .callout}

### Pull the Mandelbrot image as a SIF file on Setonix

> ## Run on Setonix
>
> Run the following commands in the terminal connected to Setonix. This location remains in effect until another location callout appears.
{: .callout}

Request an interactive allocation if you are not already working on a compute node:

```bash
$ salloc -N 1 -n 1 -c 4 --reservation=ContainersTraining -t 4:00:00
```
{: .source}

Load Pawsey's MPI-enabled Singularity module:

```bash
$ module load singularity/4.1.0-mpi
```
{: .source}

Create your personal image-library directory if needed:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```
{: .source}

Set the same Docker ID used when publishing the image and reconstruct the registry-qualified image reference in the Setonix shell:

```bash
$ DOCKER_ID="<docker-id>"
$ MPI_REMOTE_IMAGE="docker.io/${DOCKER_ID}/mandelbrot-mpi:2026.09"
```
{: .source}

Pull the Docker/OCI image and give the resulting SIF file an explicit name:

```bash
$ singularity pull \
    "${MY_LOCAL_LIBRARY}/mandelbrot-mpi--2026.09.sif" \
    "docker://${MPI_REMOTE_IMAGE}"
```
{: .source}

Singularity retrieves the manifest and filesystem layers from the registry, assembles their final filesystem state, and creates a read-only SIF image.

Define the image path and inspect the file:

```bash
$ SINGULARITY_MPI_IMAGE="${MY_LOCAL_LIBRARY}/mandelbrot-mpi--2026.09.sif"
$ ls -lh "$SINGULARITY_MPI_IMAGE"
$ singularity inspect "$SINGULARITY_MPI_IMAGE"
```
{: .source}

Inspect the preserved image-build inputs:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" \
    ls -l /opt/build-info-and-recipes/mandelbrot-mpi
```
{: .source}

Verify the packaged application without starting an MPI job:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" mpi-mandelbrot --help
```
{: .source}

Inspect the third-party acknowledgements retained in the image:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" \
    cat /usr/local/share/doc/mpi-mandelbrot/THIRD_PARTY_NOTICES.md
```
{: .source}

Now we can submit a job in Setonix that uses this image. In Setonix, move to the mandelbrot demo directory:

```bash
$ cd $TUTO/demos/build_madelbrot_docker
```
{: .source}

The Setonix script `runMandelbrotSingularityPawsey.slurm.sh` requests 16 Slurm tasks. Its defaults deliberately use a different centre and a larger workload than the local Docker test:

```text
WIDTH=6000
HEIGHT=4000
ITERATIONS=2000
CENTRE_REAL=-0.743643887037151
CENTRE_IMAGINARY=0.131825904205330
SCALE=0.002
```
{: .output}

Submit it with:

```bash
$ sbatch runMandelbrotSingularityPawsey.slurm.sh
```
{: .source}

The host-side `srun` command starts one `singularity exec` per Slurm task, following the Pawsey hybrid MPI model covered in the MPI container episode. The larger workload gives the 16 ranks substantially more pixel and iteration work than the default local test.

After the job finishes:

> ## Run on your local computer
>
> Return to a terminal on your local computer. The following commands run locally, not on Setonix.
{: .callout}

Define the Pawsey username explicitly because the local username may differ from the Pawsey username:

```bash
$ PAWSEY_USER="<pawsey-username>"
```
{: .source}

Copy the generated PNG from Setonix into the current local directory:

```bash
$ scp \
    "${PAWSEY_USER}@setonix.pawsey.org.au:/scratch/courses01/${PAWSEY_USER}/singularity-containers/demos/build_mandelbrot_docker/output/mandelbrot.singularity.setonix.png" \
    .
```
{: .source}

The file can now be opened with the local operating system's normal image viewer.

The Docker/OCI image has now been built, tested locally, published through a registry, converted into a SIF image, launched through the supported Setonix MPI model, and its result transferred back to the local computer.

### Confirm the complete workflow

Using the commands from this episode, identify the artefact or service produced at each stage:

1. `docker build --platform linux/amd64 --file mandelbrot_mpi.dockerfile --tag "$MPI_IMAGE" .`
2. `docker push "$MPI_REMOTE_IMAGE"`
3. `singularity pull OUTPUT.sif docker://REGISTRY/NAMESPACE/IMAGE:TAG`
4. `singularity run OUTPUT.sif`

> ## Solution
>
> 1. `docker build` creates a local layered Docker/OCI image in Docker's image store.
> 2. `docker push` uploads that image to a registry under the registry-qualified name and tag.
> 3. `singularity pull` retrieves the Docker/OCI image and converts it into the explicitly named SIF file.
> 4. `singularity run` starts a container process from the SIF image and executes its default action.
{: .solution}

#### Optional: transfer an image without a registry

A registry is normally the simplest and most traceable distribution method. If a registry cannot be used, Docker can export the local image to an archive.

> ## Run on your local computer
>
> Run the following archive-creation and transfer commands on the local computer where the Docker image was built.
{: .callout}

Create the Docker archive:

```bash
$ docker image save \
    --output mandelbrot-mpi--2026.09.tar \
    "$MPI_IMAGE"
```
{: .source}

Define the Pawsey username explicitly because the local username may differ from the Pawsey username:

```bash
$ PAWSEY_USER="<pawsey-username>"
```
{: .source}

Transfer the archive directly to the personal Singularity image-library directory on Setonix:

```bash
$ scp \
    mandelbrot-mpi--2026.09.tar \
    "${PAWSEY_USER}@setonix.pawsey.org.au:/software/projects/courses01/${PAWSEY_USER}/singularity/images/"
```
{: .source}

The destination directory must already exist.

> ## If the destination directory does not exist
>
> **On Setonix**, log in and create the personal image-library directory:
>
> ```bash
> $ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
> $ mkdir -p "$MY_LOCAL_LIBRARY"
> ```
> {: .source}
>
> **Back on your local computer**, run the `scp` command.
{: .solution}

After the archive has been transferred:

> ## Run on Setonix
>
> Run the remaining archive-conversion commands in the terminal connected to Setonix.
{: .callout}

Load the MPI-enabled Singularity module, define the image-library path, and move into that directory:

```bash
$ module load singularity/4.1.0-mpi
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ cd "$MY_LOCAL_LIBRARY"
```
{: .source}

Confirm that the transferred archive is present:

```bash
$ ls -lh mandelbrot-mpi--2026.09.tar
```
{: .source}

Convert the Docker archive into a SIF image:

```bash
$ singularity build \
    mandelbrot-mpi--2026.09.sif \
    docker-archive://mandelbrot-mpi--2026.09.tar
```
{: .source}

The `docker-archive://` source identifies an archive created by `docker image save`. `singularity build` reads its image layers and metadata and creates the SIF file. The Docker archive itself is not a SIF file.

Verify the resulting image and confirm that the updated Mandelbrot options are available:

```bash
$ ls -lh mandelbrot-mpi--2026.09.sif
$ singularity exec \
    mandelbrot-mpi--2026.09.sif \
    mpi-mandelbrot --help
```
{: .source}

After verifying the SIF, the transferred Docker archive can be removed:

```bash
$ rm mandelbrot-mpi--2026.09.tar
```
{: .source}

> ## Optional: load the archive into another Docker installation
>
> If the archive is transferred to another computer running Docker instead of Setonix, load it into that Docker installation with:
>
> ```bash
> $ docker image load --input mandelbrot-mpi--2026.09.tar
> ```
> {: .source}
{: .solution}

A Docker image archive may be substantially larger than a compressed SIF. Prefer publishing the image to an appropriate registry and pulling it with Singularity when a registry is available.

> ## Private and restricted images
>
> Private registry access requires authentication and site-specific handling of credentials. Software licences may also restrict whether an image can be shared, exported, or executed on another system. Follow the registry, licence, project, and Pawsey security requirements that apply to the image. Do not use a public registry merely to avoid configuring an approved private distribution method.
{: .callout}

### Manage local Docker images

> ## Run on your local computer
>
> Return to the local computer where Docker is installed. The commands in this section operate on Docker's local image store, not on Setonix.
{: .callout}

Docker stores built and pulled images in its **local image store**. This is distinct from a remote registry such as Docker Hub.

If `MPI_IMAGE` is not defined in this local shell, define it again:

```bash
$ MPI_IMAGE="mandelbrot-mpi:2026.09"
```
{: .source}

List local image references:

```bash
$ docker image ls
```
{: .source}

Inspect image metadata, including architecture, labels, and runtime configuration:

```bash
$ docker image inspect "$MPI_IMAGE"
```
{: .source}

Remove a local image reference when it is no longer needed:

```bash
$ docker image rm "$MPI_IMAGE"
```
{: .source}

Removing a local image reference does not delete the corresponding repository or tag from Docker Hub. An image may also have several local names or tags that refer to the same underlying image data. If a container still refers to an image, remove that container deliberately before removing the image rather than forcing the operation.

Remove dangling images after reviewing Docker's confirmation prompt:

```bash
$ docker image prune
```
{: .source}

Avoid `docker image prune --all` unless you understand that it can remove any image not currently used by a container. Image-layer inspection, disk-usage analysis, and build-cache cleanup are covered in the advanced Docker episode.

#### Review the Docker-to-HPC workflow

In this episode, you extended the workflow introduced in the earlier episodes:

1. Existing images can be found in registries and run with Singularity.
2. When no suitable image exists, a Dockerfile can describe a customised Docker/OCI image.
3. Docker can build and test that image away from the shared HPC system.
4. A versioned, architecture-appropriate image can be published to a registry.
5. Singularity can pull the published image into a named, read-only SIF file.
6. The SIF image can then be used with the host-directory, overlay, MPI, GPU, and scheduler workflows covered elsewhere in this training.

For research or production use, preserve the Dockerfile and build inputs, record the image digest, retain the tested SIF where appropriate, and validate the final image on the target HPC system.
