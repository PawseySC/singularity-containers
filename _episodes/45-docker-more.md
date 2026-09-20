---
title: "Engineering robust Docker images for Setonix"
teaching: 35
exercises: 30
questions:
- "How can I reduce the size and software surface of a compiled application image?"
- "Why must an image intended for Setonix work with an arbitrary non-root identity and a read-only SIF filesystem?"
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
# Use Ubuntu for both compilation and runtime
FROM docker.io/ubuntu:24.04

# Install the compiler and remove temporary package data in the same layer
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Copy and compile the application inside the final image
COPY hello.cpp /tmp/hello.cpp
RUN g++ -O2 -Wall -Wextra -Wpedantic \
        -o /usr/local/bin/hello.exe \
        /tmp/hello.cpp

# Run the application with a default message
ENTRYPOINT ["hello.exe"]
CMD ["Hello from the single-stage image"]
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
# Build stage: install the compiler and compile the application
FROM docker.io/ubuntu:24.04 AS build

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    rm -rf /var/lib/apt/lists/*

COPY hello.cpp /tmp/hello.cpp
RUN mkdir -p /out \
    && g++ -O2 -Wall -Wextra -Wpedantic \
        -o /out/hello.exe \
        /tmp/hello.cpp

# Runtime stage: start again from Ubuntu without the build tools
FROM docker.io/ubuntu:24.04

# Install only the utility required for the later inspection exercise
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends file; \
    rm -rf /var/lib/apt/lists/*

# Copy only the compiled executable from the build stage
COPY --from=build /out/hello.exe /usr/local/bin/hello.exe

# Run the application with a default message
ENTRYPOINT ["hello.exe"]
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

### Design the image's runtime interface

The previous Docker episode used `CMD` alone to define a default action that users could replace completely. The examples above use `ENTRYPOINT` together with `CMD` to make the image behave like a dedicated command-line application. `ENTRYPOINT` is not required for every image, and using it is not inherently better than using `CMD` alone.

#### Use `CMD` for a replaceable default command

With `CMD` alone, the complete default command is replaced by anything supplied after the image reference. For example:

```dockerfile
CMD ["hello.exe", "Default message"]
```
{: .source}

If an image were configured this way, running it without another command would use the complete default from `CMD`:

```bash
$ docker run --rm IMAGE
```
{: .source}

Supplying a command after the image reference would replace that complete `CMD`, rather than pass an argument to `hello.exe`:

```bash
$ docker run --rm IMAGE cat /etc/os-release
```
{: .source}

This design is appropriate when the image provides a useful default action but users should be able to replace the complete command easily. It is also useful for general software environments that contain several commands rather than one primary application.

#### Use `ENTRYPOINT` for a dedicated application image

When an image represents one primary application, `ENTRYPOINT` can define that executable and `CMD` can supply its default arguments:

```dockerfile
ENTRYPOINT ["hello.exe"]
CMD ["Hello from the multi-stage image"]
```
{: .source}

The examples in this episode therefore behave as follows:

```bash
$ docker run --rm hello-hpc:multi
$ docker run --rm hello-hpc:multi "A different message"
```
{: .source}

The first command runs `hello.exe` with the default message from `CMD`. In the second command, `"A different message"` replaces `CMD`, while `ENTRYPOINT` remains `hello.exe`. This provides an application-like interface in which values after the image reference are normally arguments to the packaged application.

The trade-off is that running an unrelated command requires the executable itself to be overridden explicitly:

```bash
$ docker run --rm \
    --entrypoint /bin/bash \
    hello-hpc:multi
```
{: .source}

Use `CMD` alone when users should be able to replace the complete default command conveniently. Use `ENTRYPOINT` with `CMD` when the image represents a particular application and values after the image reference should normally be passed to that application.

Use exec form for both instructions unless shell interpretation is required. Exec form preserves argument boundaries and gives the application a clearer process and signal model. Wrapper scripts are appropriate when runtime setup is genuinely required, but they should finish with `exec "$@"` or `exec application ...` so that the application becomes the container's main process.

The final interface must be tested with both Docker and Singularity. For HPC jobs, `singularity exec` is often clearer than relying on an image's default action because the job script records the exact executable and arguments.

### Design for non-root execution

Docker often runs a container as root by default. Singularity on Setonix normally runs the process with the invoking user's host identity, and the SIF filesystem is read-only. Therefore, a successful root-based Docker test is not sufficient.

Containers share the host kernel, so runtime privileges matter. Running an application as root gives it more authority than most applications require and can increase the consequences of an application vulnerability or configuration error, particularly when writable host directories, devices or additional privileges are exposed to the container. Running the application with a non-root identity applies the principle of least privilege. For this Pawsey workflow, it is also a portability requirement because an image that works only as root with Docker will not match normal Singularity execution on Setonix.

The application should:

- read its packaged software from immutable locations such as `/usr/local`;
- write results to the current directory or an explicitly bind-mounted output directory;
- use `/tmp` or another writable location for temporary files;
- avoid writing caches or configuration under `/usr`, `/opt`, or a root-owned home directory;
- avoid `sudo`, system-user creation, service startup, and privileged initialisation at runtime;
- tolerate a runtime numeric user ID that is not listed in the image's `/etc/passwd` file.

> ## Why is the Docker daemon a separate security concern?
>
> Traditional Docker installations use a daemon that normally runs with host root privileges. The daemon performs host-level operations such as starting containers, mounting filesystems, configuring networks and exposing devices. Access to the daemon or its socket must therefore be treated as highly privileged host access.
>
> This daemon privilege is separate from the identity used by an application inside a particular container. Running an application as a non-root user reduces that application's authority, but it does not make unrestricted access to the Docker daemon safe. This distinction is one reason Docker is used for building on a developer-controlled computer while Singularity is used as the user-facing runtime on Setonix.
{: .solution}

#### Optional: define a non-root user for Docker

Images intended primarily for Docker or Kubernetes sometimes define a known non-root user so that the application does not run as root by default:

```dockerfile
RUN groupadd --system application \
    && useradd --system \
        --gid application \
        --create-home \
        application

USER application
```
{: .source}

The `USER` instruction sets the default user for subsequent `RUN` instructions and for container execution. The application and its runtime directories must be accessible to that user.

Pawsey application images do not generally need to define a fixed runtime user. Singularity normally runs the application with the invoking user's host UID and GID rather than changing to the user declared by `USER`. For Setonix portability, the application should instead tolerate an arbitrary non-root numeric identity that may not have an entry in the image's `/etc/passwd` file.

#### Test an arbitrary non-root identity with Docker

On macOS or Linux, run the image with the current host UID and primary GID. Bind mount the host working directory at `/work`, then use `hello.exe` inside the container to create a file there:

```bash
$ docker run --rm \
    --user "$(id -u):$(id -g)" \
    --mount type=bind,source="$PWD",target=/work \
    --workdir /work \
    --entrypoint /bin/bash \
    hello-hpc:multi \
    -c 'hello.exe "Created by the container" > container-output.txt'
```
{: .source}

The `--user` option runs the container process with the host user's numeric UID and primary GID. On a native Linux Docker host, this normally prevents files written through the bind mount from being owned by root or by a fixed user defined in the image. The host directory's normal permissions still apply. Access that depends on supplementary groups may require additional group configuration.

The `--entrypoint /bin/bash` option temporarily replaces the image's normal `hello.exe` entrypoint. Bash interprets the output redirection, while `hello.exe` produces the content. Because `/work` is the bind-mounted host working directory, `container-output.txt` is created in `$PWD` on the host rather than in the image.

On Windows PowerShell, use a known numeric UID and GID for the Linux-container test:

```powershell
PS> docker run --rm --user "1000:1000" --mount "type=bind,source=${PWD},target=/work" --workdir /work --entrypoint /bin/bash hello-hpc:multi -c 'hello.exe "Created by the container" > container-output.txt'
```
{: .source}

Inspect the file from the host and confirm that the host user can modify and remove it:

```bash
$ ls -ln container-output.txt
$ cat container-output.txt
$ printf '%s\n' "Modified by the host" >> container-output.txt
$ rm container-output.txt
```
{: .source}

This checks both non-root execution and usable ownership and permissions for output written to the host. Docker Desktop on macOS and Windows shares files through a Linux virtual machine, so its UID and GID behaviour may not exactly match a native Linux host or Setonix. Final validation must still be performed on Setonix.

#### Compare with Singularity on Setonix

The Docker command above tests locally how the application is expected to behave when the image is later converted to SIF and run with Singularity. The equivalent command on Setonix would be:

```bash
$ singularity exec \
    hello-hpc--multi.sif \
    bash -c 'hello.exe "Created by the container" > container-output.txt'
```
{: .source}

Singularity normally runs the command with the invoking user's host identity and makes the host current working directory available inside the container at the same path. Therefore, this example does not require Docker's `--user`, `--mount`, or `--workdir` options. In both commands, Bash interprets the redirection, `hello.exe` generates the content, and `container-output.txt` is stored in the host working directory rather than in the container image.

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

### Inspect executable and library compatibility

Check the executable inside the image:

```bash
$ docker run --rm --entrypoint /usr/bin/file hello-hpc:multi /usr/local/bin/hello.exe
$ docker run --rm --entrypoint /usr/bin/ldd hello-hpc:multi /usr/local/bin/hello.exe
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
