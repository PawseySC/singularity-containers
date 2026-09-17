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
- "Read and write a basic Dockerfile using `FROM`, `RUN`, `ENV`, `WORKDIR`, `USER`, and `CMD`"
- "Build and test an `amd64` Docker/OCI image"
- "Apply basic practices for build contexts, package installation, image tags, and runtime users"
- "Publish an image to Docker Hub and pull it as a named SIF file on an HPC system"
keypoints:
- "Docker is commonly used to build and test Docker/OCI images on a workstation, while Singularity runs the resulting images on the HPC system"
- "A Dockerfile records the base image and the instructions used to assemble a new image"
- "Docker image layers support build caching, but Dockerfile instruction order and cleanup affect build efficiency and image size"
- "Use a small build context, a `.dockerignore` file, a trusted base image, and a meaningful image tag"
- "Build for the CPU architecture of the target HPC system"
- "Do not store passwords, access tokens, or other secrets in a Dockerfile, build argument, or image layer"
- "A registry provides the normal bridge between a Docker build environment and Singularity on an HPC system"
- "Singularity can pull a Docker/OCI image from a registry and convert it into a named SIF file"
---

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

Singularity is the container engine used to run images throughout this training because it is designed for shared HPC systems and integrates with host filesystems, schedulers, MPI libraries, GPUs, and high-speed interconnects.

Docker is widely used to build Docker/OCI images on personal computers, workstations, cloud systems, and CI/CD services. Its layered image format and build cache make repeated development builds convenient. Docker/OCI images can also be distributed through standard registries and consumed by several container engines, including Singularity.

The tools therefore have complementary roles in this training:

- **Docker builds and tests the image** on a system where Docker is available.
- **A container registry distributes the image** between systems.
- **Singularity converts and runs the image** on the HPC system.

Docker is not used to run the workload on Setonix. The final execution still uses Singularity and follows the cluster-specific practices introduced in the other episodes.

> ## Build privileges and available alternatives
>
> A conventional Docker Engine installation uses a privileged daemon, although Docker Desktop and rootless Docker provide different deployment models. In all cases, build images only on systems where container building is supported and authorised.
>
> Shared HPC login and compute nodes normally do not provide Docker to users. Singularity can also build images from definition files, including with supported remote-build or fakeroot configurations, but those workflows are separate from the Docker-based workflow in this episode.
{: .callout}

#### Prepare the Docker build directory

Perform the Docker sections of this episode on the computer where Docker was installed and tested in the installation episode. Do not run these commands on a Setonix login or compute node.

Move to the Docker example in the training repository:

```bash
$ cd "$TUTO/demos/lolcow_docker"
$ pwd
```
{: .source}

The directory should contain a `Dockerfile`:

```bash
$ ls -la
```
{: .source}

A Docker build uses a **build context**. The context is the directory, URL, or other source made available to the builder. In this episode, the final `.` in the build command selects the current directory as the build context.

Keep the context small. Files in it can be sent to the builder and may become available to `COPY` and `ADD` instructions. Unrelated data, Git history, generated output, credentials, and large files should not be included.

Create a `.dockerignore` file for common files that the example does not need:

```text
.git
.gitignore
*.sif
*.tar
*.tar.gz
```
{: .source}

A `.dockerignore` file works similarly to `.gitignore`: matching paths are excluded from the build context.

#### Read the Dockerfile

The Dockerfile for the example is:

```dockerfile
FROM ubuntu:24.04

LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre"

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        fortune-mod \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/games:${PATH}"

RUN useradd --create-home --uid 1000 training
USER training
WORKDIR /home/training

CMD ["bash", "-c", "fortune | cowsay | lolcat"]
```
{: .source}

Docker reads the instructions from top to bottom. Each instruction describes part of the resulting image or its default runtime configuration.

##### `FROM`: select a base image

```dockerfile
FROM ubuntu:24.04
```
{: .source}

`FROM` begins a build stage and selects its base image. This example starts from the versioned `ubuntu:24.04` image rather than `ubuntu:latest`.

Choose base images from trusted publishers and prefer a supported, suitably small image that provides what the application needs. A versioned tag communicates the intended base more clearly, although tags can still be updated by their publisher. For stricter provenance, production builds may pin the base image by digest and update that digest deliberately.

##### `LABEL`: record image metadata

```dockerfile
LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre"
```
{: .source}

`LABEL` adds metadata to the image. OCI annotation names are used here so that the purpose and publisher of the image can be identified with image-inspection tools.

##### `RUN`: execute build-time commands

```dockerfile
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        fortune-mod \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*
```
{: .source}

`RUN` executes commands while the image is being built and stores the resulting filesystem changes in an image layer.

This instruction:

1. enables strict and verbose shell behaviour for the build step
2. updates the Ubuntu package index
3. installs the three required packages without additional recommended packages
4. removes cached package data that is not needed at runtime

The package-index update, installation, and cleanup are performed in the same `RUN` instruction. Removing files in a later layer would not remove them from the earlier layer in which they were created.

Using fewer `RUN` instructions does not by itself guarantee a good Dockerfile. Combine commands when their changes belong in one filesystem layer, but keep the result readable and ensure failures stop the build.

##### `ENV`: define a runtime environment variable

```dockerfile
ENV PATH="/usr/games:${PATH}"
```
{: .source}

`ENV` defines an environment variable that persists in the image and is normally present when a container is started. Ubuntu installs `cowsay` and `fortune` under `/usr/games`, so this instruction adds that directory to `PATH`.

An `export` performed inside one `RUN` instruction affects only the shell used for that build step. Use `ENV` when the setting should form part of the image's runtime environment.

##### `USER`: avoid running the application as root

```dockerfile
RUN useradd --create-home --uid 1000 training
USER training
```
{: .source}

Build steps that install operating-system packages require root inside the Docker build. The image then creates an unprivileged account and uses `USER` to select it for subsequent instructions and for Docker containers started from the image.

This is a useful default for testing an image with Docker. On the HPC system, Singularity applies its own runtime identity model and normally runs the containerised process as the invoking cluster user.

##### `WORKDIR`: select the default directory

```dockerfile
WORKDIR /home/training
```
{: .source}

`WORKDIR` sets the working directory for later Dockerfile instructions and for the default container execution. Unlike `RUN cd ...`, it remains in effect after the build step has ended.

The working directory selected by a Dockerfile does not prevent Singularity from starting in a host directory that it bind mounts into the container. As seen in the basic Singularity episode, Singularity normally starts a containerised command in the host current working directory.

##### `CMD`: define the default action

```dockerfile
CMD ["bash", "-c", "fortune | cowsay | lolcat"]
```
{: .source}

`CMD` defines the default command used when a Docker container is started without another command. The JSON, or exec, form preserves the command and arguments as separate values.

The pipeline must be interpreted by a shell, so the Dockerfile explicitly starts `bash -c`. This is the same shell-expression principle used with `singularity exec` in the basic Singularity episode.

> ## `CMD` is a default, not a build command
>
> `RUN` executes while the image is built. `CMD` records the default command to execute later when a container is started. A Dockerfile can contain only one effective `CMD`; if several are present, the last one takes effect.
{: .callout}

#### Build the image for the target HPC architecture

Set a meaningful image tag for the exercise:

```bash
$ IMAGE_TAG="lolcow:2026.09"
```
{: .source}

Setonix compute nodes use the `amd64` architecture, also called `x86_64`. Build explicitly for that target:

```bash
$ docker build --platform linux/amd64 -t "$IMAGE_TAG" .
```
{: .source}

The `-t` option assigns the repository name `lolcow` and tag `2026.09`. The final `.` selects the current directory as the build context.

On an `amd64` Linux or Windows system, the requested architecture matches the host. On an Apple Silicon or other `arm64` system, Docker Desktop normally uses emulation to perform this build, so it may take longer.

> ## If the image is not loaded into the local Docker image store
>
> Some Docker installations use a `buildx` builder whose output is not loaded automatically. In that case, build with:
>
> ```bash
> $ docker buildx build --platform linux/amd64 --load -t "$IMAGE_TAG" .
> ```
> {: .source}
{: .callout}

List the local image:

```bash
$ docker image ls "$IMAGE_TAG"
```
{: .source}

Inspect the architecture recorded for the image:

```bash
$ docker image inspect "$IMAGE_TAG" --format '{{.Os}}/{{.Architecture}}'
```
{: .source}

The expected output is:

```text
linux/amd64
```
{: .output}

### Build the image again

Run the same build command a second time:

```bash
$ docker build --platform linux/amd64 -t "$IMAGE_TAG" .
```
{: .source}

Which steps are reused from the build cache? Why is the second build normally faster?

> ## Solution
>
> Docker can reuse layers when the instruction and the files on which it depends have not changed. Because neither the Dockerfile nor the build context changed, most or all build steps should report that cached results were used.
>
> If an early instruction changes, Docker must rebuild that step and later steps that depend on it. Dockerfile instruction order therefore affects how effectively the build cache can be reused.
{: .solution}

#### Test the image with Docker

Run the image's default command:

```bash
$ docker run --rm "$IMAGE_TAG"
```
{: .source}

The output should contain a fortune displayed by `cowsay` and coloured by `lolcat`. The exact message and colours vary between runs.

The `--rm` option removes the stopped container after it exits. It does not remove the image.

Override the default `CMD` by supplying another command after the image name:

```bash
$ docker run --rm "$IMAGE_TAG" id
```
{: .source}

The output should identify the unprivileged `training` user rather than root.

Open an interactive shell for inspection:

```bash
$ docker run --rm -it "$IMAGE_TAG" bash
```
{: .source}

The options have the following roles:

- `-i` keeps standard input open
- `-t` allocates a terminal
- `--rm` removes the stopped container when the shell exits

Inside the container, inspect the environment:

```bash
training@CONTAINER-ID:~$ cat /etc/os-release
training@CONTAINER-ID:~$ command -v fortune cowsay lolcat
training@CONTAINER-ID:~$ pwd
training@CONTAINER-ID:~$ exit
```
{: .source}

The container's writable layer is temporary in this example because `--rm` removes the container after it stops. Changes made interactively do not update the image or the Dockerfile. If an interactive test reveals a required change, edit the Dockerfile and rebuild the image.

This is the reproducible development loop:

1. edit the Dockerfile or required build-context files
2. build the image
3. test the image
4. repeat until the Dockerfile reliably produces the required environment

#### Basic practices for research and HPC images

A working Dockerfile is only the starting point. When preparing an image for a research workflow:

- Start from an image maintained by a trusted project, vendor, or organisation.
- Prefer explicit application and base-image versions over `latest` where practical.
- Record the Dockerfile and related build files in version control.
- Keep the build context small and use `.dockerignore`.
- Install only required packages and remove package-manager caches in the same layer.
- Put frequently changing instructions and copied source files after stable dependency-installation steps when practical, so the build cache can be reused.
- Use `COPY` for ordinary file and directory copies. Use `ADD` only when its additional behaviour is specifically required.
- Use an unprivileged runtime user unless the application has a documented reason not to.
- Do not place passwords, private keys, access tokens, licence files, or other secrets in the Dockerfile, build arguments, environment variables, or copied build context. Use the build system's supported secret mechanism when a build must access protected resources.
- Rebuild and test images regularly so that base-image and package security updates are incorporated.
- Test the final SIF image on the target HPC system, including MPI, GPU, filesystem, and performance behaviour where relevant.

For compiled applications, consider a **multi-stage build**. A build stage can contain compilers and development packages, while a later runtime stage receives only the installed application and required runtime libraries. This can substantially reduce the final image size and software surface. The MPI episode includes an example of the build-time requirements for an HPC application linked against MPI.

#### Publish the image to Docker Hub

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
$ docker tag "$IMAGE_TAG" "$REMOTE_IMAGE"
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

#### Pull the image as a SIF file on Setonix

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

Notice that the identity differs from the Docker test. Under Docker, the Dockerfile's `USER training` instruction selected the runtime user. Under the normal Singularity execution model on the cluster, the process runs with your cluster user identity rather than becoming the image's `training` user.

### Confirm the complete workflow

Using the commands from this episode, identify the artefact or service produced at each stage:

1. `docker build --platform linux/amd64 -t "$IMAGE_TAG" .`
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
$ docker image save -o lolcow--2026.09.tar "$IMAGE_TAG"
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
