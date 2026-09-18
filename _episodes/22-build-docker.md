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
- "Publish an image to Docker Hub and pull it as a named SIF file on an HPC system"
keypoints:
- "Docker is commonly used to build and test Docker/OCI images on a workstation, while Singularity runs the resulting images on the HPC system"
- "A Dockerfile records the base image and the instructions used to assemble a new image"
- "`COPY` adds files from the build context to the image"
- "Docker image layers support build caching, but Dockerfile instruction order and cleanup affect build efficiency and image size"
- "Use a small build context, a `.dockerignore` file, a trusted base image, and a meaningful image tag"
- "Build for the CPU architecture of the target HPC system"
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
$ cd demos/lolcow_docker
$ pwd
```
{: .source}

The working directory should end with:

```text
singularity-containers/demos/lolcow_docker
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
$ cd "$TUTO/demos/lolcow_docker"
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

#### `CMD`: define the default action

```dockerfile
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]
```
{: .source}

`CMD` defines the default command used when a Docker container is started without another command. The JSON, or exec, form preserves the command and arguments as separate values.

The redirection and pipeline must be interpreted by a shell, so the Dockerfile explicitly starts `bash -c`. `cowsay` reads the copied message through standard input, and `lolcat` colours the resulting output. This is the same shell-expression principle used with `singularity exec` in the basic Singularity episode.

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

The final `.` selects the current directory as the **build context**. The build context is the collection of files and directories made available to the Docker builder. Files in the context can be used by instructions such as `COPY` and `ADD`, so the context should contain only what the build requires. Unrelated data, Git history, generated output, credentials, and large files should not be included.

Create a `.dockerignore` file to exclude common files that this build does not need:

```text
.git
.gitignore
*.sif
*.tar
*.tar.gz
```
{: .source}

A `.dockerignore` file works similarly to `.gitignore`: matching paths are excluded from the build context before it is sent to the builder.

The same recipe can be selected explicitly with the `-f` option instead of relying on the `Dockerfile` symbolic link:

```bash
$ docker build \
    --platform linux/amd64 \
    -f lolcow.dockerfile \
    -t "$IMAGE" \
    .
```
{: .source}

Both commands build from `lolcow.dockerfile`. The first demonstrates Docker's conventional default filename, while the second names the recipe explicitly.

On an `amd64` Linux or Windows system, the requested architecture matches the host. On an Apple Silicon or other `arm64` system, Docker Desktop normally uses emulation to perform this build, so it may take longer.

> ## If the image is not loaded into the local Docker image store
>
> Some Docker installations use a `buildx` builder whose output is not loaded automatically. In that case, build with:
>
> ```bash
> $ docker buildx build --platform linux/amd64 --load -t "$IMAGE" .
> ```
> {: .source}
{: .callout}

List the local image:

```bash
$ docker image ls "$IMAGE"
```
{: .source}

Inspect the architecture recorded for the image:

```bash
$ docker image inspect "$IMAGE" --format '{{.Os}}/{{.Architecture}}'
```
{: .source}

The expected output is:

```text
linux/amd64
```
{: .output}

### Modify the copied message and rebuild the image

Display the message file on the local computer:

```bash
$ cat lolcow-message.txt
```
{: .source}

Its initial content is:

```text
Built with Docker and ready to run with Singularity!
```
{: .output}

Replace the message with one of your own:

```bash
$ printf '%s\n' 'My updated container message' > lolcow-message.txt
```
{: .source}

Rebuild the image using the same name and tag:

```bash
$ docker build --platform linux/amd64 -t "$IMAGE" .
```
{: .source}

Which build steps are reused from the cache, and which step runs again?

> ## Solution
>
> Docker can reuse the unchanged base-image and package-installation layers. The `COPY` step runs again because `lolcow-message.txt`, one of its inputs, changed. Instructions after that changed step are also reconsidered.
>
> Placing stable and expensive installation steps before frequently changing application files allows Docker to reuse more of the build cache during development.
{: .solution}

Run the rebuilt image and confirm that it displays the new message:

```bash
$ docker run --rm "$IMAGE"
```
{: .source}

#### Test the image with Docker

Run the image's default command:

```bash
$ docker run --rm "$IMAGE"
```
{: .source}

The output should contain the message from `lolcow-message.txt`, displayed by `cowsay` and coloured by `lolcat`.

The `--rm` option removes the stopped container after it exits. It does not remove the image.

Override the default `CMD` by supplying another command after the image name:

```bash
$ docker run --rm "$IMAGE" cat /etc/os-release
```
{: .source}

Docker runs the selected command instead of the image's default `CMD`.

Open an interactive shell for inspection:

```bash
$ docker run --rm -it "$IMAGE" bash
```
{: .source}

The options have the following roles:

- `-i` keeps standard input open
- `-t` allocates a terminal
- `--rm` removes the stopped container when the shell exits

Inside the container, inspect the environment:

```bash
root@CONTAINER-ID:/# cat /etc/os-release
root@CONTAINER-ID:/# command -v cowsay lolcat
root@CONTAINER-ID:/# pwd
root@CONTAINER-ID:/# exit
```
{: .source}

The container's writable layer is temporary in this example because `--rm` removes the container after it stops. Changes made interactively do not update the image or the Dockerfile. If an interactive test reveals a required change, edit the Dockerfile and rebuild the image.

This is the reproducible development loop:

1. edit the Dockerfile or required build-context files
2. build the image
3. test the image
4. repeat until the Dockerfile reliably produces the required environment

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
