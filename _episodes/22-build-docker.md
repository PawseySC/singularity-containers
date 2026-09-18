---
title: "Building Docker/OCI images for HPC with Singularity"
teaching: 30
exercises: 20
questions:
- "Why are Dockerfiles commonly used to build images that will run with Singularity on HPC systems?"
- "How can I describe and build a reproducible container image with a Dockerfile?"
- "How can I test the image locally before using it on a cluster?"
- "How can I publish an image to a registry and convert it to a SIF file?"
objectives:
- "Explain the roles of Docker and Singularity in an HPC container workflow"
- "Read and write a basic Dockerfile using `FROM`, `LABEL`, `RUN`, `ENV`, `COPY`, and `CMD`"
- "Build and test an `amd64` Docker/OCI image"
- "Apply basic practices for build contexts, package installation, image tags, and runtime users"
- "Build an MPI application image from a Pawsey-provided base image"
- "Publish an image to Docker Hub and pull it as a named SIF file on an HPC system"
keypoints:
- "Docker is commonly used to build and test Docker/OCI images on a workstation, while Singularity runs the resulting images on the HPC system"
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

Perform the Docker sections of this episode on the local computer where Docker was installed and tested in the installation episode. Do not run these Docker commands on a Setonix login or compute node.

Open a terminal on your local computer and clone the training repository if you have not already done so:

```bash
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
$ cd "$TUTO"
```
{: .source}

Now move to the working directory for this episode:

```bash
$ cd demos/build_lolcow_docker
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

The tools therefore have complementary roles in this training:

- **Docker builds and tests the Docker/OCI image** on the participant's local computer.
- **A container registry stores and distributes the image** in a widely supported format.
- **Singularity retrieves, converts, and runs the image** on the HPC system.

Docker is not used to run the workload on Setonix. The final execution uses Singularity and follows the cluster-specific practices introduced in the other episodes.

### Prepare the Docker build directory

Perform the Docker sections of this episode on the computer where Docker was installed and tested in the installation episode. Do not run these commands on a Setonix login or compute node.

Move to the Docker example in the training repository:

```bash
$ cd "$TUTO/demos/build_lolcow_docker"
$ pwd
```
{: .source}

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

The `lolcow-message.txt` file contains the message that will be copied into the image and displayed by the default container action.

### Read the Dockerfile

The `lolcow.dockerfile` recipe is:

```dockerfile
# Start from a versioned Ubuntu image on Docker Hub
FROM docker.io/ubuntu:24.04

# Record standard OCI image metadata
LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre"

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

#### `LABEL`: record image metadata

```dockerfile
LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre"
```
{: .source}

`LABEL` adds metadata to the image and accepts one or more key-value pairs. This example uses predefined annotation keys from the [OCI Image Specification](https://specs.opencontainers.org/image-spec/annotations/), including `org.opencontainers.image.title`, `org.opencontainers.image.description`, and `org.opencontainers.image.vendor`. Using these standard keys makes the metadata easier for OCI-compatible tools to interpret consistently.

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

#### `CMD`: define the default action

```dockerfile
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]
```
{: .source}

`CMD` defines the default command used when a Docker container is started without another command. The JSON, or exec, form preserves the command and arguments as separate values.

The redirection and pipeline must be interpreted by a shell, so the Dockerfile explicitly starts `bash -c`. `cowsay` reads the copied message through standard input, and `lolcat` processes the resulting output. This is the same shell-expression principle used with `singularity exec` in the basic Singularity episode.

Later in the episode, the `CMD` instruction is updated to use `lolcat --force`, ensuring that colour codes are emitted even when Docker has not allocated a terminal.

> ## `CMD` is a default, not a build command
>
> `RUN` executes while the image is built. `CMD` records the default command to execute later when a container is started. A Dockerfile can contain only one effective `CMD`; if several are present, the last one takes effect.
{: .callout}

### Build the image for the target HPC architecture

Define the complete local image reference, including its repository name and tag:

```bash
$ IMAGE="lolcow:2026.09"
```
{: .source}

Setonix compute nodes use the `amd64` architecture, also called `x86_64`. Build explicitly for that target:

```bash
$ docker build --platform linux/amd64 -t "$IMAGE" .
```
{: .source}

The `IMAGE` variable contains the local repository name `lolcow` and tag `2026.09`. The `-t` option assigns that complete reference to the image. Docker finds the recipe through the symbolic link named `Dockerfile`.

The final `.` selects the current directory as the **build context**. The build context is the collection of files and directories made available to the Docker builder. Files in the context can be used by instructions such as `COPY` and `ADD`, so the context should contain only what the build requires.

> ## Limiting a larger build context
>
> This example has a small, controlled build context containing only the recipe and its message file, so it does not require a `.dockerignore` file. Larger projects commonly use `.dockerignore` to reduce unnecessary build-context processing and to prevent broad `COPY` or `ADD` instructions from including unwanted or sensitive files.
{: .callout}

List the newly built image:

```bash
$ docker image ls "$IMAGE"
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
$ docker image inspect "$IMAGE" --format '{% raw %}{{.Os}}/{{.Architecture}}{% endraw %}'
```
{: .source}

The expected output is:

```text
linux/amd64
```
{: .output}

Now run the image's default action without `--rm`:

```bash
$ docker run "$IMAGE"
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
$ docker run --rm "$IMAGE"
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
> $ docker run --rm --platform linux/amd64 "$IMAGE"
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
$ docker build --platform linux/amd64 --tag "$IMAGE" .
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
$ docker build \
    --platform linux/amd64 \
    --file lolcow.dockerfile \
    --tag "$IMAGE" \
    .
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
$ docker buildx build \
    --check \
    --file lolcow.dockerfile \
    .
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
$ docker build --platform linux/amd64 -t "$IMAGE" .
```
{: .source}

Which build steps are reused from the cache, and which part changes?

> ## Solution
>
> The base image, package-installation layer, environment setting, and copied message are unchanged, so Docker can reuse their cached results. Only the final image configuration changes because the `CMD` instruction was modified.
>
> This demonstrates why instructions that change frequently are normally placed after stable and expensive build steps.
{: .solution}

Run the rebuilt image:

```bash
$ docker run --rm "$IMAGE"
```
{: .source}

The same message should now be displayed in colour. The terminal must support ANSI colour sequences, as standard terminals on current Linux, macOS, and Windows installations normally do.

### Run other commands with Docker

Override the default `CMD` by supplying another command after the image name:

```bash
$ docker run --rm "$IMAGE" cat /etc/os-release
```
{: .source}

Docker runs the selected command instead of the image's default `CMD`.

### Open an interactive shell with Docker

As with `singularity shell`, Docker can open an interactive shell for inspecting and testing an image:

```bash
$ docker run --rm --interactive --tty "$IMAGE" bash
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
    "$IMAGE" \
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
    "$IMAGE" \
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

### Basic practices for research and HPC images

A working Dockerfile is only the starting point. When preparing an image for a research workflow:

- Start from an image maintained by a trusted project, vendor, or organisation.
- Prefer explicit application and base-image versions over `latest` where practical.
- Record the Dockerfile and related build files in version control.
- Keep the build context small and use `.dockerignore`.
- Install only required packages and remove package-manager caches in the same layer.
- Put frequently changing `COPY` instructions after stable dependency-installation steps when practical, so the build cache can be reused.
- Use `COPY` for ordinary file and directory copies. Use `ADD` only when its additional behaviour is specifically required.
- Do not place passwords, private keys, access tokens, licence files, or other secrets in the Dockerfile, build arguments, environment variables, or copied build context. Use the build system's supported secret mechanism when a build must access protected resources.
- Rebuild and test images regularly so that base-image and package security updates are incorporated.
- Test the final SIF image on the target HPC system, including MPI, GPU, filesystem, and performance behaviour where relevant.

For compiled applications, consider a **multi-stage build**. A build stage can contain compilers and development packages, while a later runtime stage receives only the installed application and required runtime libraries. This can substantially reduce the final image size and software surface. The MPI episode includes an example of the build-time requirements for an HPC application linked against MPI.

### Use Pawsey-provided base images

Pawsey publishes container base images that users can extend for their own applications. These provide tested starting environments for Pawsey systems, including MPI base images designed for compatibility with the Cray MPICH environment on Setonix and ROCm-based images prepared for AMD GPU workloads. A derived application image and its complete workflow must still be tested on the target system.

The recipes used to build Pawsey-supported images are available from the [Pawsey container recipes repository](https://github.com/PawseySC/pawsey-containers). The corresponding Docker/OCI images are published under the [Pawsey organisation on Quay.io](https://quay.io/pawsey).

For this example, move to a separate build context:

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
render-mandelbrot
```
{: .output}

We will build an MPI Mandelbrot renderer on top of this Pawsey-provided image:

```text
quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04
```
{: .output}

The base image already provides MPICH, the GNU compiler toolchain, and MPI utilities. The derived image compiles an MPI C++ application with the base image's `mpic++` compiler and installs ImageMagick to convert the application's raw PPM result into a PNG file that is straightforward to view.

The training application was independently written for this example and informed by the educational MPI partitioning approaches in Liam Ryan's MIT-licensed Mandelbrot repository. Attribution and the upstream licence are retained in `mpi-mandelbrot.cpp` and `THIRD_PARTY_NOTICES.md`. The notice file is also copied into the final image.

#### Read the MPI application recipe

The complete `mandelbrot_mpi.dockerfile` recipe is:

```dockerfile
# Build the application on Pawsey's Setonix-compatible MPICH base image
FROM quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04

# Record standard OCI image metadata
LABEL org.opencontainers.image.title="MPI Mandelbrot renderer" \
      org.opencontainers.image.description="MPI training application built on Pawsey's MPICH base image" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre" \
      org.opencontainers.image.licenses="MIT"

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

# Install the wrapper that launches MPI and converts the result to PNG
COPY render-mandelbrot /usr/local/bin/render-mandelbrot
RUN chmod 0755 /usr/local/bin/render-mandelbrot

# Preserve third-party acknowledgements and licence information
COPY THIRD_PARTY_NOTICES.md \
    /usr/local/share/doc/mpi-mandelbrot/THIRD_PARTY_NOTICES.md

# Display the wrapper help when no other command is supplied
CMD ["render-mandelbrot", "--help"]
```
{: .source}

This recipe demonstrates how a specialised base image can provide build tools and runtime libraries. The application is compiled with `mpic++` from the Pawsey base image, so the derived image uses the MPI environment supplied and tested by Pawsey.

The `render-mandelbrot` wrapper launches the MPI program, writes a temporary PPM image, converts it to PNG, and removes the temporary file. The C++ source divides image rows among MPI ranks and uses `MPI_Gatherv` to assemble the calculated pixels on rank 0.

#### Build the MPI application image

Define its local image reference:

```bash
$ MPI_IMAGE="mandelbrot-mpi:2026.09"
```
{: .source}

Build it for Setonix's CPU architecture:

```bash
$ docker build \
    --platform linux/amd64 \
    --file mandelbrot_mpi.dockerfile \
    --tag "$MPI_IMAGE" \
    .
```
{: .source}

Confirm that the image exists:

```bash
$ docker image ls "$MPI_IMAGE"
```
{: .source}

Run the default action to display the wrapper help:

```bash
$ docker run --rm --platform linux/amd64 "$MPI_IMAGE"
```
{: .source}

#### Render a Mandelbrot image with MPI

Run four MPI processes and bind mount the current host directory at `/work` so that the PNG result persists after the container exits:

```bash
$ docker run \
    --rm \
    --platform linux/amd64 \
    --mount type=bind,source="$PWD",target=/work \
    "$MPI_IMAGE" \
    render-mandelbrot \
        --processes 4 \
        --width 1200 \
        --height 800 \
        --iterations 500 \
        --output /work/mandelbrot.png
```
{: .source}

Representative output is:

```text
MPI Mandelbrot renderer
Image size: 1200 x 800
Maximum iterations: 500
MPI processes: 4
PPM output: /tmp/tmp.XXXXXXXXXX.ppm
Rendering completed in 0.420 seconds
PNG output: /work/mandelbrot.png
```
{: .output}

The temporary filename and elapsed time will differ. Confirm that the PNG file exists on the host:

```bash
$ ls -lh mandelbrot.png
```
{: .source}

Open `mandelbrot.png` with the normal image viewer or web browser on the local computer.

The MPI execution in this section occurs entirely within Docker on one local computer. Running the resulting MPI image across Setonix compute nodes requires the Singularity, Slurm, and host-MPI integration covered in the MPI container episode.

### Publish the image to Docker Hub

A registry is the normal way to move a Docker/OCI image from the build computer to the HPC system. This section uses Docker Hub, which was introduced in the basic Singularity episode.

You need a Docker Hub account and a repository to push the image. Replace `<dockerhub-account>` with your account name:

```bash
$ DOCKERHUB_ACCOUNT="<dockerhub-account>"
$ REMOTE_IMAGE="docker.io/${DOCKERHUB_ACCOUNT}/lolcow:2026.09"
```
{: .source}

Authenticate from the Docker client:

```bash
$ docker login docker.io
```
{: .source}

Follow the authentication instructions shown by Docker. Do not put the password or access token directly in the Dockerfile or shell history.

Add the registry-qualified tag to the existing local image:

```bash
$ docker tag "$IMAGE" "$REMOTE_IMAGE"
```
{: .source}

The new tag is another name for the same local image. Confirm both names:

```bash
$ docker image ls
```
{: .source}

Push the registry-qualified image:

```bash
$ docker push "$REMOTE_IMAGE"
```
{: .source}

After the push completes, inspect the repository and tag in Docker Hub. For this training exercise, the repository must be readable from the cluster. Do not publish proprietary, confidential, export-controlled, licensed, or otherwise restricted software in a public repository.

> ## Image tags and immutable identity
>
> A tag such as `2026.09` is a readable reference, but a registry publisher can move a tag to different image content. A pushed image is also identified by a content digest beginning with `sha256:`.
>
> For a reproducible workflow, record the registry, namespace, repository, tag, digest, Dockerfile revision, and relevant build inputs. A retained SIF file also preserves the exact content that was pulled at that time.
{: .callout}

### Pull the image as a SIF file on Setonix

Return to the terminal connected to Setonix and request an interactive allocation if you are not already working on a compute node:

```bash
$ salloc -N 1 -n 1 -c 4 --reservation=ContainersTraining -t 4:00:00
```
{: .source}

Load the Singularity module used for non-MPI containers:

```bash
$ module load singularity/4.1.0-nohost
```
{: .source}

Create your personal image-library directory if needed:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```
{: .source}

Set the same Docker Hub account name used when publishing the image:

```bash
$ DOCKERHUB_ACCOUNT="<dockerhub-account>"
```
{: .source}

Pull the Docker/OCI image and give the resulting SIF file an explicit name:

```bash
$ singularity pull \
    "${MY_LOCAL_LIBRARY}/lolcow--2026.09.sif" \
    "docker://docker.io/${DOCKERHUB_ACCOUNT}/lolcow:2026.09"
```
{: .source}

Singularity retrieves the manifest and filesystem layers from the registry, assembles their final filesystem state, and creates a read-only SIF image.

Define the image path and inspect the file:

```bash
$ SINGULARITY_IMAGE="${MY_LOCAL_LIBRARY}/lolcow--2026.09.sif"
$ ls -lh "$SINGULARITY_IMAGE"
$ singularity inspect "$SINGULARITY_IMAGE"
```
{: .source}

Run the default action recorded by the Docker image:

```bash
$ singularity run "$SINGULARITY_IMAGE"
```
{: .source}

Select another command with `singularity exec`:

```bash
$ singularity exec "$SINGULARITY_IMAGE" id
```
{: .source}

Under the normal Singularity execution model on the cluster, the containerised process runs with your cluster user identity. The Docker image does not need to contain a user account matching each possible HPC user.

### Confirm the complete workflow

Using the commands from this episode, identify the artefact or service produced at each stage:

1. `docker build --platform linux/amd64 -t "$IMAGE" .`
2. `docker push "$REMOTE_IMAGE"`
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

A registry is normally the simplest and most traceable distribution method. If a registry cannot be used, Docker can export a local image to an archive:

```bash
$ docker image save -o lolcow--2026.09.tar "$IMAGE"
```
{: .source}

The archive can be transferred to another Docker installation and loaded with:

```bash
$ docker image load -i lolcow--2026.09.tar
```
{: .source}

A Docker image archive is not a SIF file. It preserves the Docker/OCI image representation for Docker-compatible tooling and may be substantially larger than a compressed SIF. For the workflow taught here, prefer publishing the image to an appropriate registry and pulling it with Singularity.

> ## Private and restricted images
>
> Private registry access requires authentication and site-specific handling of credentials. Software licences may also restrict whether an image can be shared, exported, or executed on another system. Follow the registry, licence, project, and Pawsey security requirements that apply to the image. Do not use a public registry merely to avoid configuring an approved private distribution method.
{: .callout}

#### Review the Docker-to-HPC workflow

In this episode, you extended the workflow introduced in the earlier episodes:

1. Existing images can be found in registries and run with Singularity.
2. When no suitable image exists, a Dockerfile can describe a customised Docker/OCI image.
3. Docker can build and test that image away from the shared HPC system.
4. A versioned, architecture-appropriate image can be published to a registry.
5. Singularity can pull the published image into a named, read-only SIF file.
6. The SIF image can then be used with the host-directory, overlay, MPI, GPU, and scheduler workflows covered elsewhere in this training.

For research or production use, preserve the Dockerfile and build inputs, record the image digest, retain the tested SIF where appropriate, and validate the final image on the target HPC system.
