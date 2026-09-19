---
title: "Engineering robust Docker images for Setonix"
teaching: 35
exercises: 30
questions:
- "How can I reduce the size and software surface of a compiled application image?"
- "How can I design an image that works without root privileges and with a read-only SIF filesystem?"
- "How should `ENTRYPOINT`, `CMD`, data, configuration, and launch policy be separated?"
- "How can I inspect binaries, image layers, disk usage, and build cache before validating an image on Setonix?"
objectives:
- "Create and explain a Docker multi-stage build"
- "Compare image size and contents before and after optimisation"
- "Test an image with a non-root user and writable bind-mounted paths"
- "Design a portable runtime interface using `ENTRYPOINT`, `CMD`, arguments, and environment variables"
- "Inspect executable architecture and shared-library dependencies"
- "Manage local images and Docker build cache safely"
- "Apply the checks to MPI and GPU images intended for Setonix"
keypoints:
- "A Docker multi-stage build uses multiple `FROM` instructions and `COPY --from` to separate build tools from the runtime image"
- "A smaller image comes from excluding unnecessary content, not simply from minimising the number of layers"
- "Setonix containers must work with the user's host identity and writable host directories"
- "Images should contain applications and stable dependencies, while data, credentials, configuration, and Slurm launch policy remain outside"
- "`file`, `ldd`, image inspection, and target-system testing reveal compatibility problems that a successful Docker build cannot"
- "Images, containers, build cache, and remote registry content are separate objects and must be managed deliberately"
---

### Scope of this episode

This episode follows [Building Docker/OCI images for HPC with Singularity]({% link _episodes/22-build-docker.md %}). That episode built, tested, published, and ran an image on Setonix. Here, the goal is to turn a working Dockerfile into a smaller, safer, and more predictable image for research and HPC use.

> ## Run on your local computer
>
> Run the Docker commands in this episode on the local computer where Docker is installed. Do not run Docker builds on a Setonix login or compute node. A later section identifies the validation commands that run on Setonix.
{: .callout}

Set the training repository path and create a separate working directory:

```bash
$ export TUTO="$PWD/singularity-containers"
$ mkdir -p "$TUTO/demos/docker_advanced"
$ cd "$TUTO/demos/docker_advanced"
```
{: .source}

### Reduce image size and software surface

The first topic in this episode is reducing the size and software surface of a compiled application image. We will first build a single-stage image that contains both the build tools and the application, then convert it to a multi-stage build whose final stage contains only the application and its required runtime dependencies. Comparing the two images will show what a multi-stage build removes and what must remain for the application to run correctly.

#### Build a single-stage compiled image

Create `hello.cpp`:

```cpp
#include <iostream>

int main(int argc, char **argv) {
    const char *message = argc > 1 ? argv[1] : "Hello from a Setonix-ready image";
    std::cout << message << '\n';
    return 0;
}
```
{: .source}

Create `hello.single.dockerfile`:

```dockerfile
FROM docker.io/ubuntu:24.04

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

COPY hello.cpp /tmp/hello.cpp
RUN g++ -O2 -Wall -Wextra -Wpedantic \
        -o /usr/local/bin/hello \
        /tmp/hello.cpp

CMD ["hello"]
```
{: .source}

Build it for Setonix's CPU architecture:

```bash
$ docker build \
    --platform linux/amd64 \
    --file hello.single.dockerfile \
    --tag hello-hpc:single \
    .
```
{: .source}

Run and inspect it:

```bash
$ docker run --rm --platform linux/amd64 hello-hpc:single
$ docker image ls hello-hpc:single
$ docker image history hello-hpc:single
```
{: .source}

The image works, but it also contains the compiler, development files, source file, and package-management content that remain after the build.

#### Convert the recipe to a multi-stage build

A Docker multi-stage build contains multiple `FROM` instructions. `AS build` gives the first stage a stable name, and `COPY --from=build` copies selected artefacts into the final stage.

Create `hello.multistage.dockerfile`:

```dockerfile
FROM docker.io/ubuntu:24.04 AS build

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    rm -rf /var/lib/apt/lists/*

COPY hello.cpp /tmp/hello.cpp
RUN mkdir -p /out \
    && g++ -O2 -Wall -Wextra -Wpedantic \
        -o /out/hello \
        /tmp/hello.cpp

FROM docker.io/ubuntu:24.04

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends file; \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /out/hello /usr/local/bin/hello

ENTRYPOINT ["hello"]
CMD ["Hello from the multi-stage image"]
```
{: .source}

Build and run the final stage:

```bash
$ docker build \
    --platform linux/amd64 \
    --file hello.multistage.dockerfile \
    --tag hello-hpc:multi \
    .
$ docker run --rm --platform linux/amd64 hello-hpc:multi
```
{: .source}

Compare the size and layer history of the two final images:

```bash
$ docker image ls hello-hpc
$ docker image history hello-hpc:single
$ docker image history hello-hpc:multi
```
{: .source}

The output from `docker image ls` should show that `hello-hpc:multi` is smaller than `hello-hpc:single`. The history of the single-stage image includes the installation of `g++` and the compilation steps because they form part of that final image. In contrast, the history of the multi-stage runtime image begins from its second `FROM` instruction and includes only the runtime-stage instructions. The build stage is not included in the final runtime image.

The multi-stage image therefore excludes the compiler, development files, and source code while retaining the compiled executable and its required runtime libraries. `docker image history` shows how the final image was constructed; it is not a complete inventory of the files or packages inside the image.

> ## MPI and GPU runtime stages
>
> Do not copy an MPI or GPU executable into an arbitrary small base image merely to reduce size. The final stage must provide ABI-compatible runtime libraries, and the image must remain compatible with the Pawsey Singularity module and host-library injection used by the workload. Start from a Pawsey-provided base image where applicable and validate the result on Setonix.
{: .callout}

#### Keep the build context small with `.dockerignore`

Create `.dockerignore`:

```text
.git
*.tar
*.sif
output/
```
{: .source}

A `.dockerignore` file excludes matching paths from the build context sent to Docker. This reduces build-context processing and helps prevent unrelated or sensitive files from being available to broad `COPY` or `ADD` instructions. It does not remove content that a Dockerfile explicitly creates or downloads during the build, and it reduces the final image size only when the excluded content would otherwise have been copied into the image.

#### Understand what actually reduces image size

Useful size-reduction practices include:

1. Choose an appropriate supported base image.
2. Install only packages required at runtime.
3. Use `--no-install-recommends` where appropriate.
4. Remove package indexes and temporary build files in the same `RUN` instruction that creates them.
5. Keep the build context small and exclude unnecessary files with `.dockerignore`.
6. Use multi-stage builds to leave compilers, headers, source files, and intermediate artefacts out of the final stage.

As discussed in the previous Docker episode, files created in one image layer remain part of that layer even if a later instruction deletes them. For operations such as package installation, keep the package-index update, installation, and cleanup in the same `RUN` instruction so that temporary package data is not retained in an earlier layer.

This does not mean that every command should be concatenated into a single `RUN` instruction. Combine commands whose filesystem changes belong together, while keeping separate logical build steps readable and allowing Docker to reuse useful cached layers. Fewer layers do not automatically produce a smaller or better image.

### Design for non-root execution

Docker often runs a container as root by default. Singularity on Setonix normally runs the process with your host user identity, and the SIF filesystem is read-only. Therefore, a successful root-based Docker test is not sufficient.

The application should:

- read its packaged software from immutable locations such as `/usr/local`;
- write results to the current directory or an explicitly bind-mounted output directory;
- use `/tmp` or another writable location for temporary files;
- avoid writing caches or configuration under `/usr`, `/opt`, or a root-owned home directory;
- avoid `sudo`, system-user creation, service startup, and privileged initialisation at runtime;
- tolerate a runtime numeric user ID that is not listed in the image's `/etc/passwd` file.

Test the image with the current host UID and GID on macOS or Linux:

```bash
$ docker run --rm \
    --user "$(id -u):$(id -g)" \
    --mount type=bind,source="$PWD",target=/work \
    --workdir /work \
    hello-hpc:multi \
    "Running as an unprivileged user"
```
{: .source}

On Windows PowerShell, use a known numeric UID and GID for a Linux-container test:

```powershell
PS> docker run --rm --user "1000:1000" --mount "type=bind,source=${PWD},target=/work" --workdir /work hello-hpc:multi "Running as an unprivileged user"
```
{: .source}

A stronger application test should also create its expected output in `/work` and confirm that the host user can read, modify, and remove that file afterward.

### Separate the image from runtime data and launch policy

Use this design model:

```text
Image                 application and stable runtime dependencies
Bind mounts           input data and writable output directories
Arguments and ENV     runtime configuration
Host script           local execution policy
Slurm job script      Setonix resources and launch policy
Registry or SIF       image distribution and immutable execution artefact
```
{: .output}

Do not bake user data, project paths, credentials, Slurm reservations, process counts, or site-specific launch commands into the image.

### Design `ENTRYPOINT` and `CMD` deliberately

In exec form:

- `ENTRYPOINT` defines the executable that should normally run.
- `CMD` supplies default arguments, or supplies the default command when no `ENTRYPOINT` is defined.
- arguments after the image reference in `docker run` replace `CMD` while retaining `ENTRYPOINT`.

For the multi-stage example:

```bash
$ docker run --rm hello-hpc:multi
$ docker run --rm hello-hpc:multi "A different message"
```
{: .source}

Avoid shell form unless shell interpretation is required. Exec form preserves argument boundaries and gives the application a clearer process and signal model. Wrapper scripts are appropriate when runtime setup is genuinely required, but they should finish with `exec "$@"` or `exec application ...` so that the application becomes the container's main process.

The final interface must be tested with both Docker and Singularity. For HPC jobs, `singularity exec` is often clearer than relying on an image's default action because the job script records the exact executable and arguments.

### Inspect executable and library compatibility

Check the executable inside the image:

```bash
$ docker run --rm --entrypoint /usr/bin/file hello-hpc:multi /usr/local/bin/hello
$ docker run --rm --entrypoint /usr/bin/ldd hello-hpc:multi /usr/local/bin/hello
```
{: .source}

`file` should report an x86-64 executable for a `linux/amd64` build. `ldd` lists the shared libraries that must be present in the runtime stage. A line containing `not found` indicates a missing runtime dependency.

For an MPI image, also record the MPI implementation and inspect the application dependencies:

```bash
$ docker run --rm "$MPI_IMAGE" mpiexec --version
$ docker run --rm "$MPI_IMAGE" ldd /usr/local/bin/mpi-mandelbrot
```
{: .source}

A successful local build cannot prove compatibility with Setonix MPI, GPU, scheduler, or injected host libraries. Those features must be validated with the appropriate Pawsey base image, Singularity module, and Slurm launch model on Setonix.

### Manage images, containers, disk usage, and build cache

List local image references and include digests where available:

```bash
$ docker image ls
$ docker image ls --digests
```
{: .source}

Inspect an image and its layer history:

```bash
$ docker image inspect hello-hpc:multi
$ docker image history hello-hpc:multi
```
{: .source}

Inspect Docker disk usage:

```bash
$ docker system df
```
{: .source}

Images, stopped containers, and build cache are separate objects. Remove them deliberately:

```bash
$ docker container ls --all
$ docker container prune
$ docker image prune
$ docker builder prune
```
{: .source}

Each prune command asks for confirmation. `docker image prune` removes dangling images by default, while `docker builder prune` removes unused build cache and may make later builds slower. Avoid broad cleanup commands such as `docker system prune --all` in shared or valuable development environments unless you have first reviewed exactly what they will remove.

Remove a specific local reference with:

```bash
$ docker image rm hello-hpc:single
```
{: .source}

Removing a local reference does not delete a tag from Docker Hub. Conversely, deleting a remote repository does not clean the local image store.

### Build secrets and private dependencies

Do not pass passwords, tokens, private keys, or licence files with `ARG`, `ENV`, or ordinary `COPY`. Those mechanisms can preserve sensitive material in image metadata, layers, or build cache.

When a build must access a protected resource, use BuildKit secret or SSH mounts. For example:

```bash
$ docker build \
    --secret id=service_token,src="$HOME/.service-token" \
    --file protected.dockerfile \
    .
```
{: .source}

The corresponding Dockerfile instruction makes the secret available only to that build step:

```dockerfile
RUN --mount=type=secret,id=service_token \
    application-that-reads /run/secrets/service_token
```
{: .source}

Secrets must still be excluded from the build context and `.dockerignore` should be treated as a secondary safeguard, not as the secret-delivery mechanism.

### Validate the final image on Setonix

> ## Run on Setonix
>
> Publish the final `linux/amd64` image or transfer it through the approved workflow, pull or convert it to a named SIF file, and run the following checks in an interactive allocation using the Singularity module appropriate for the workload.
{: .callout}

At minimum, verify:

1. `singularity inspect` reports the expected metadata.
2. The packaged executable starts with `singularity exec`.
3. Inputs and outputs use writable Setonix paths.
4. The process runs with your Setonix identity and does not require root.
5. `file` and `ldd` show the expected architecture and available runtime libraries.
6. MPI and GPU applications use the supported Pawsey base image, Singularity module, and Slurm launch model.
7. A representative job produces correct output before performance measurements are trusted.

### Review the robust-image workflow

A robust Setonix image is not merely one that completes `docker build`. It should be small enough to distribute efficiently, contain only required runtime software, work without root privileges, keep data and site policy outside the image, expose a predictable runtime interface, and pass validation with Singularity on the target system.
