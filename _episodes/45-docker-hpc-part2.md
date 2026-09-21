---
title: "Building Docker/OCI images for HPC, Part 2: Robust and portable images"
teaching: 35
exercises: 30
questions:
- "How can I reduce the size of an image?"
- "What practices make a Docker image robust and suitable for HPC use?"
- "How can I validate a Docker image on my computer before running it on Setonix?"
- "How can private dependencies be accessed during a build without storing credentials in the image?"
- "How can I manage local Docker images, containers and build cache, and distinguish them from images stored in a registry?"
objectives:
- "Create and explain a Docker multi-stage build"
- "Compare the size, construction history and contents of single-stage and multi-stage images"
- "Keep the Docker build context small with `.dockerignore`"
- "Choose an appropriate runtime interface using `CMD` and `ENTRYPOINT`"
- "Test an image with an arbitrary non-root identity and writable bind-mounted paths"
- "Separate applications and stable dependencies from data, configuration and launch policy"
- "Inspect executable architecture and shared-library dependencies with Docker and Singularity"
- "Explain how Pawsey's MPI-enabled Singularity environment changes runtime library resolution"
- "Access private repositories and packages using BuildKit SSH and secret mounts"
- "Manage local images, containers and build cache safely"
- "Validate representative MPI and GPU workflows on Setonix"
keypoints:
- "Multi-stage builds leave compilers, source files and other build-only content out of the final runtime image"
- "Files created in one layer remain in that layer even if a later instruction deletes them"
- "A `.dockerignore` file controls the build context but does not remove content created by Dockerfile instructions"
- "`CMD` provides a replaceable default command, while `ENTRYPOINT` can define the executable of a dedicated application image"
- "Images intended for Setonix should tolerate an arbitrary non-root identity and write persistent data through host directories"
- "Applications and stable dependencies belong in the image, while data, credentials, runtime configuration and site-specific launch policy remain outside"
- "BuildKit SSH and secret mounts provide temporary access to protected resources without copying credentials into image layers"
- "`file` and `ldd` can reveal architecture and dependency problems but cannot prove MPI, GPU or application correctness"
- "Pawsey's MPI-enabled Singularity module makes compatible Cray MPICH and Setonix communication libraries available to containerised MPI applications"
- "A representative Slurm job on Setonix remains the final correctness and integration test"
- "Docker images, stopped containers, build cache, registry content and SIF files are separate objects that must be managed deliberately"
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

The image works, but it also contains the compiler, development files, source file and package-management content that remain after the build.

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

The multi-stage image therefore excludes the compiler, development files and source code while retaining the compiled executable and its required runtime libraries. `docker image history` shows how the final image was constructed; it is not a complete inventory of the files or packages inside the image.

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
6. Use multi-stage builds to leave compilers, headers, source files and intermediate artefacts out of the final stage.

As discussed in the previous Docker episode, files created in one image layer remain part of that layer even if a later instruction deletes them. For operations such as package installation, keep the package-index update, installation and cleanup in the same `RUN` instruction so that temporary package data is not retained in an earlier layer.

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
- avoid `sudo`, system-user creation, service startup and privileged initialisation at runtime;
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

A container image should contain the software environment that remains stable between executions. It should not permanently contain information that changes for each dataset, user, system or job.

For example, the Mandelbrot image from the previous Docker episode contains:

- the `mpi-mandelbrot` executable;
- the MPI implementation and runtime libraries;
- ImageMagick, which converts the generated PPM file to PNG;
- licence information and image metadata.

These components belong in the image because every execution requires them and they should travel together as one tested software environment.

Other parts of the workflow remain outside the image because they depend on the particular execution:

- input data and output directories;
- image dimensions, iteration count, centre and scale;
- the number of MPI processes;
- local filesystem paths;
- Slurm resources, reservations and time limits;
- whether Docker, `mpiexec`, Singularity or `srun` launches the application;
- registry credentials and other secrets.

Use the following separation as a design guide:

```text
Image                 application and stable runtime dependencies
Bind mounts           input data and writable output directories
Arguments and ENV     runtime configuration
Host script           local execution policy
Slurm job script      Setonix resources and launch policy
Registry or SIF       image distribution and immutable execution artefact
```
{: .output}

#### Application and stable dependencies belong in the image

The image should contain the executable and the libraries and utilities required wherever that application runs. Keeping these components together provides the consistent software environment that the container is intended to preserve.

Installing dependencies or compiling the application at runtime would make each execution slower and less predictable. These operations should normally occur while the image is built.

#### Input and output belong on the host

Research data should normally remain outside the image. Docker and Singularity bind mounts make host directories available to the container without copying their contents permanently into the image.

This allows the same image to process different datasets and write results to writable project or scratch directories. Rebuilding the image should not be necessary merely because an input file, output location or dataset has changed.

#### Arguments and environment variables configure a run

Values that change between executions should normally be supplied at runtime. Command-line arguments are appropriate for explicit application options, while environment variables are useful for configuration that the application or host-side launch script reads from its environment.

For the Mandelbrot example, values such as these configure a particular rendering:

```bash
WIDTH=3000
HEIGHT=2000
ITERATIONS=1000
CENTRE_REAL=-0.743643887037151
CENTRE_IMAGINARY=0.131825904205330
SCALE=0.002
```
{: .source}

These values should not be fixed permanently in the image because changing the rendering should not require rebuilding and republishing the software environment.

#### Host scripts describe local launch policy

A host-side script records how an image should be used in a particular environment. For a local Docker test, the script can select the image, bind mount an output directory, choose the number of local MPI processes, run `mpiexec` inside one container and convert the result to PNG.

These operations describe how the application is launched locally. They are not intrinsic properties of the application itself.

#### Slurm scripts describe Setonix launch policy

A Slurm job script contains decisions that belong to the HPC system, including:

- partition, reservation and time limit;
- node, task and CPU requests;
- the Singularity module to load;
- the path to the SIF image;
- the `srun` and `singularity exec` commands;
- project and scratch paths.

These details should not be placed in the image. If they were embedded there, changing the requested resources, reservation or launch method would require rebuilding and redistributing the image. Keeping the Slurm policy in a host-side script also makes the exact resources and launch command visible when the job is reviewed.

The same application image can therefore participate in different workflows:

```text
Local Docker:
host script -> docker run -> container mpiexec -> application ranks

Setonix:
Slurm job -> srun -> singularity exec -> application ranks
```
{: .output}

The image remains the same, but each environment provides its appropriate launcher, writable paths and resource policy. Do not bake user data, project paths, credentials, Slurm reservations, process counts or site-specific launch commands into the image.

### Inspect executable and library compatibility

A successful image build does not prove that a compiled executable has the correct architecture or that all its runtime libraries are available. The following checks use the Mandelbrot image from the previous Docker episode to inspect the packaged executable locally and then compare its library resolution on Setonix.

#### Inspect the Mandelbrot image with Docker

> ## Run on your local computer
>
> Run the following Docker commands on the local computer where the Mandelbrot image was built. Do not run them on Setonix. This location remains in effect until the **Run on Setonix** callout below.
{: .callout}

Define the image reference:

```bash
$ MPI_IMAGE="mandelbrot-mpi:2026.09"
```
{: .source}

Confirm the image architecture:

```bash
$ docker image inspect "$MPI_IMAGE" --format '{% raw %}{{.Os}}/{{.Architecture}}{% endraw %}'
```
{: .source}

The expected result is:

```text
linux/amd64
```
{: .output}

Inspect the executable format and shared-library dependencies:

```bash
$ docker run --rm \
    --platform linux/amd64 \
    --entrypoint /usr/bin/file \
    "$MPI_IMAGE" \
    /usr/local/bin/mpi-mandelbrot

$ docker run --rm \
    --platform linux/amd64 \
    --entrypoint /usr/bin/ldd \
    "$MPI_IMAGE" \
    /usr/local/bin/mpi-mandelbrot
```
{: .source}

Representative output from `file` includes:

```text
/usr/local/bin/mpi-mandelbrot: ELF 64-bit LSB pie executable, x86-64, dynamically linked, ...
```
{: .output}

The important points are that the executable is a Linux ELF binary, is built for `x86-64`, and is dynamically linked. The `ldd` output should resolve every library and must not contain `not found`. Locally, the MPI dependency resolves to the MPICH library packaged in the image:

```text
libmpi.so.12 => /usr/lib/libmpi.so.12 (...)
libstdc++.so.6 => /lib/x86_64-linux-gnu/libstdc++.so.6 (...)
libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (...)
```
{: .output}

Record the MPI implementation packaged in the image:

```bash
$ docker run --rm \
    --platform linux/amd64 \
    --entrypoint mpiexec \
    "$MPI_IMAGE" \
    --version
```
{: .source}

The Pawsey MPICH base image used in the previous episode reports MPICH 3.4.3. This documents the MPI implementation used to build the application, but does not by itself prove compatibility with Setonix.

#### Inspect the SIF image on Setonix

> ## Run on Setonix
>
> Run the following commands in an interactive allocation on a Setonix compute node. Load Pawsey's MPI-enabled Singularity module before inspecting the runtime library resolution.
{: .callout}

```bash
$ module load singularity/4.1.0-mpi
$ SINGULARITY_MPI_IMAGE="${MYSOFTWARE}/singularity/images/mandelbrot-mpi--2026.09.sif"
```
{: .source}

To inspect the executable format, read the packaged binary inside the container and pipe it to the host `file` command:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" \
    cat /usr/local/bin/mpi-mandelbrot \
    | /usr/bin/file -
```
{: .source}

The host shell interprets the pipe. `cat` runs inside the container and reads the packaged executable, while `/usr/bin/file` runs on Setonix and inspects the bytes received through standard input. This avoids mixing the container's `file` utility with host libraries injected by the MPI-enabled module.

Representative output is:

```text
/dev/stdin: ELF 64-bit LSB shared object, x86-64, dynamically linked, ...
```
{: .output}

Depending on the `file` version and its recognition database, a position-independent executable may be described as either a PIE executable or a shared object. The important results are `x86-64` and `dynamically linked`.

The Setonix MPI environment injects and preloads numerous host libraries, so its complete `ldd` output is substantially longer than the local Docker output. First, filter the output to show the libraries most relevant to MPI and interconnect integration:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" \
    ldd /usr/local/bin/mpi-mandelbrot \
    | grep -Ei 'mpi|fabric|cxi|xpmem|pmi|pals'
```
{: .source}

Representative lines include:

```text
/opt/xpmem/lib64/libxpmem.so.0 (...)
/usr/lib64/libcxi.so.1 (...)
libmpi.so.12 => /opt/cray/pe/mpich/8.1.32/ofi/gnu/12.3/lib-abi-mpich/libmpi.so.12 (...)
libfabric.so.1 => /opt/cray/libfabric/1.22.0/lib64/libfabric.so.1 (...)
libpmi.so.0 => /opt/cray/pe/pmi/default/lib/libpmi.so.0 (...)
libpmi2.so.0 => /opt/cray/pe/pmi/default/lib/libpmi2.so.0 (...)
libpals.so.0 => /opt/cray/pe/pals/default/lib/libpals.so.0 (...)
```
{: .output}

Check separately for unresolved dependencies:

```bash
$ singularity exec "$SINGULARITY_MPI_IMAGE" \
    ldd /usr/local/bin/mpi-mandelbrot \
    | grep 'not found'
```
{: .source}

No output is the expected result.

The exact versions and paths may change when the Setonix software environment is updated. The important comparison is that local Docker execution resolves `libmpi.so.12` to the generic MPICH library packaged in the image, while the MPI-enabled Singularity module redirects that compatible dependency to Setonix's Cray MPICH library and makes the host communication and process-management libraries available.

> ## Optional: inspect a larger MPI application
>
> The Mandelbrot image demonstrated the main compatibility checks. Apply the same checks to the OpenFOAM image used in the MPI episode.
>
> Confirm that:
>
> 1. `pimpleFoam` is an `x86-64`, dynamically linked ELF executable.
> 2. No shared-library dependency is reported as `not found`.
> 3. `libmpi.so.12` resolves to Setonix's Cray MPICH compatibility library.
> 4. The Setonix communication and process-management libraries are available.
> 5. The packaged OpenFOAM application starts successfully.
>
> Assume that `SINGULARITY_IMAGE` contains the path to the OpenFOAM SIF image:
>
> ```bash
> $ SINGULARITY_IMAGE="${MYSOFTWARE}/singularity/images/openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif"
> ```
> {: .source}
>
> > ## Solution
> >
> > Locate `pimpleFoam`:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" \
> >     bash -c 'command -v pimpleFoam'
> > ```
> > {: .source}
> >
> > Read the packaged executable inside the container and inspect it with the host `file` command:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" \
> >     bash -c 'cat "$(command -v pimpleFoam)"' \
> >     | /usr/bin/file -
> > ```
> > {: .source}
> >
> > Representative output includes:
> >
> > ```text
> > /dev/stdin: ELF 64-bit LSB shared object, x86-64, dynamically linked, ...
> > ```
> > {: .output}
> >
> > Inspect the most relevant runtime-library resolution:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" \
> >     bash -c 'ldd "$(command -v pimpleFoam)"' \
> >     | grep -Ei 'mpi|fabric|cxi|xpmem|pmi|pals|not found'
> > ```
> > {: .source}
> >
> > Representative lines include:
> >
> > ```text
> > /opt/xpmem/lib64/libxpmem.so.0 (...)
> > /usr/lib64/libcxi.so.1 (...)
> > libmpi.so.12 => /opt/cray/pe/mpich/8.1.32/ofi/gnu/12.3/lib-abi-mpich/libmpi.so.12 (...)
> > libfabric.so.1 => /opt/cray/libfabric/1.22.0/lib64/libfabric.so.1 (...)
> > libpmi.so.0 => /opt/cray/pe/pmi/default/lib/libpmi.so.0 (...)
> > libpmi2.so.0 => /opt/cray/pe/pmi/default/lib/libpmi2.so.0 (...)
> > libpals.so.0 => /opt/cray/pe/pals/default/lib/libpals.so.0 (...)
> > ```
> > {: .output}
> >
> > Check explicitly for unresolved dependencies:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" \
> >     bash -c 'ldd "$(command -v pimpleFoam)"' \
> >     | grep 'not found'
> > ```
> > {: .source}
> >
> > No output is the expected result.
> >
> > Finally, confirm that the application starts:
> >
> > ```bash
> > $ singularity exec "$SINGULARITY_IMAGE" \
> >     pimpleFoam -help \
> >     | head -n 30
> > ```
> > {: .source}
> >
> > The expected output identifies OpenFOAM v2606 and displays the `pimpleFoam` options.
> >
> > These checks show that the executable has the correct architecture, its libraries resolve, and Pawsey's MPI-enabled Singularity module exposes the Cray MPICH and Setonix communication libraries. They do not replace the OpenFOAM Slurm exercise. A representative parallel job using `srun`, `singularity exec` and `pimpleFoam -parallel` remains the final functional validation.
> {: .solution}
{: .challenge}

A line containing `not found` indicates an unresolved runtime dependency and must be investigated. However, successful `file` and `ldd` checks still do not prove that MPI communication works correctly. Complete the validation by running a representative Slurm job through the supported launch pattern:

```text
srun -> singularity exec -> mpi-mandelbrot
```
{: .output}

For MPI and GPU images, architecture and library inspection are starting points. Correctness, multi-node communication, hardware integration and performance must still be validated with the appropriate Pawsey base image, Singularity module and Slurm launch model on Setonix.

### Handle secrets and private dependencies securely

An image build may need authentication to retrieve private source code, packages or other protected resources. The downloaded source or installed software may belong in the image, but the credential used to obtain it must not.

Do not provide passwords, access tokens, SSH private keys or licence files through `ARG`, `ENV` or ordinary `COPY`. Sensitive information supplied through these mechanisms may remain in image metadata, filesystem layers, build history or build cache.

BuildKit provides two related mechanisms:

- SSH mounts provide temporary access to an SSH agent, commonly for cloning private Git repositories.
- Secret mounts provide a token, password, licence file or other credential temporarily as a file or environment variable.

#### Clone a private Git repository with an SSH mount

Do not copy an SSH private key into the image:

```dockerfile
COPY id_ed25519 /root/.ssh/id_ed25519
RUN git clone git@example.org:research/private-code.git
```
{: .source}

Even if a later instruction removes the key, it remains stored in the layer created by `COPY`.

Instead, keep the private key on the local computer and make it available through the local SSH agent. Pass access to that agent into the build:

```bash
$ docker build \
    --ssh default \
    --file private-code.dockerfile \
    .
```
{: .source}

The Dockerfile can request the SSH mount only for the instruction that clones the repository:

```dockerfile
RUN --mount=type=ssh \
    git clone \
        git@example.org:research/private-code.git \
        /tmp/private-code
```
{: .source}

The build process can authenticate through the forwarded SSH agent, but the private key itself remains on the local computer and is not copied into the resulting image layer.

In a multi-stage build, the private repository would normally be cloned and compiled in the build stage. Only the required executable or other runtime artefacts would then be copied into the final stage.

#### Access a private package registry with a secret mount

A private package registry or API may use an access token instead of SSH. Keep the token in a protected file on the local computer, for example:

```text
$HOME/.service-token
```
{: .output}

Pass the token separately from the ordinary build context:

```bash
$ docker build \
    --secret id=service_token,src="$HOME/.service-token" \
    --file private-package.dockerfile \
    .
```
{: .source}

The Dockerfile can temporarily mount the token for the instruction that retrieves the protected package:

```dockerfile
RUN --mount=type=secret,id=service_token \
    curl \
        --fail \
        --header "Authorization: Bearer $(cat /run/secrets/service_token)" \
        --output /tmp/private-package.tar.gz \
        https://packages.example.org/private-package.tar.gz \
    && tar -xzf /tmp/private-package.tar.gz -C /opt \
    && rm /tmp/private-package.tar.gz
```
{: .source}

During this instruction, the token is available at:

```text
/run/secrets/service_token
```
{: .output}

The secret mount disappears when the `RUN` instruction finishes. The installed package remains, but the secret mount is not stored in the resulting layer. The command must still avoid copying the secret into the image or printing it in the build output.

#### Keep credentials out of the build context

A `.dockerignore` file should exclude private keys, token files and other credentials from the ordinary build context:

```text
id_ed25519
id_rsa
.service-token
*.key
*.pem
```
{: .source}

The mechanisms serve different purposes:

```text
.dockerignore             prevents credential files from entering
                          the ordinary build context

--ssh                     temporarily forwards SSH-agent access

--secret                  temporarily mounts a token, password,
                          licence file or another credential
```
{: .output}

A file excluded by `.dockerignore` can still be supplied through `--secret`, because the secret mount uses a separate BuildKit channel rather than the ordinary build context.

Before placing private source code, licensed software or protected packages in an image, confirm that the applicable licence and project policies permit the image to be stored, transferred or published. Never push the resulting image to a public registry unless that distribution is explicitly permitted.

### Validate the final image on Setonix

> ## Run on Setonix
>
> Final validation must use the Singularity module and launch model appropriate for the workload. Perform these checks in an interactive allocation or through a representative Slurm job on a Setonix compute node.
{: .callout}

Before treating an image as ready for research workloads, verify that:

1. `singularity inspect` reports the expected image metadata.
2. The packaged executable has the expected `x86-64` architecture.
3. All required shared libraries resolve, with no dependencies reported as `not found`.
4. Applications intended for MPI or GPU use resolve the supported Pawsey host libraries under the appropriate Singularity module.
5. The application runs with the invoking user's Setonix identity and does not require root.
6. Inputs, outputs, temporary files and caches use writable host locations rather than the read-only SIF filesystem.
7. Runtime configuration is supplied through arguments, environment variables or host-side configuration rather than being fixed unnecessarily in the image.
8. Slurm resource requests and site-specific launch policy remain in the job script.
9. A representative job produces correct output using the supported Slurm and Singularity launch pattern.
10. Performance is measured only after correctness, library integration and multi-node behaviour have been established.

The `file` and `ldd` checks identify obvious architecture and dependency problems, but they do not prove that MPI communication, GPU access or application behaviour is correct. The final validation is a representative workload launched through the supported Setonix execution model.

### Manage local Docker storage and registry references

> ## Run on your local computer
>
> Return to the local computer where Docker is installed. The commands in this section operate on Docker's local storage, not on Setonix or on a remote registry.
{: .callout}

Docker manages several related but separate objects:

- **Images** are the read-only content used to create containers. One image can have several local names or tags.
- **Containers** are runnable instances created from images. A stopped container can remain after its main process exits unless it was started with `--rm`.
- **Build cache** stores intermediate build results that can accelerate later builds.
- **Registry repositories and tags** are remote references stored in services such as Docker Hub or Quay.io. They are not part of the local Docker image store.

#### Review local storage before removing anything

List local image references and stopped or running containers:

```bash
$ docker image ls
$ docker container ls --all
```
{: .source}

Inspect Docker's summary of local disk usage:

```bash
$ docker system df
```
{: .source}

Use `docker image inspect IMAGE` when you need detailed metadata for a particular image, and `docker image history IMAGE` when you need to review its recorded construction steps. These commands inspect the local image and do not contact or modify a registry.

#### Remove specific objects when possible

Prefer removing a known container or image reference rather than beginning with a broad prune command:

```bash
$ docker container rm CONTAINER
$ docker image rm hello-hpc:single
```
{: .source}

An image can have more than one local tag. Removing one tag removes that reference, but Docker retains the underlying image data while another tag or container still refers to it. If a container still uses the image, remove the container deliberately before removing the final image reference rather than forcing the operation.

A locally removed image can normally be recovered by rebuilding it from its Dockerfile or pulling it again from a registry, provided the required source and registry content remain available.

#### Prune unused objects deliberately

After reviewing the local state, these commands provide narrower cleanup operations:

```bash
$ docker container prune
$ docker image prune
$ docker builder prune
```
{: .source}

Each command asks for confirmation:

- `docker container prune` removes stopped containers.
- `docker image prune` removes dangling image data by default.
- `docker builder prune` removes unused build cache and may make subsequent builds slower.

Avoid broad cleanup commands such as `docker system prune --all` unless you have reviewed exactly what they can remove. In routine work, specific removals and the narrower prune commands are easier to reason about.

#### Local and registry content are independent

A tag such as:

```text
docker.io/example-user/mandelbrot-mpi:2026.09
```
{: .output}

can identify both a local image reference and a remote registry tag, but those are stored separately. The following operations affect only the local Docker installation:

```bash
$ docker image rm docker.io/example-user/mandelbrot-mpi:2026.09
$ docker builder prune
```
{: .source}

They do not delete the remote tag from Docker Hub or another registry. Conversely, deleting a repository or tag through a registry does not remove local images, stopped containers, build cache or SIF files already created from that image.

Deleting remote content requires the registry's own interface or API and should be performed only after confirming that other users and workflows no longer require it. Keep important Dockerfiles, build inputs, version tags and image digests recorded so that local and remote artefacts can be identified and reproduced.

### Review the robust-image workflow

A robust image is not merely one that completes `docker build`. The complete workflow is:

1. Keep the build context focused and exclude unnecessary or sensitive files.
2. Use multi-stage builds to separate compilation from the final runtime environment.
3. Retain every runtime library required by the packaged executable.
4. Choose `CMD` alone or `ENTRYPOINT` with `CMD` according to the intended runtime interface.
5. Ensure that the application works with an arbitrary non-root identity and a read-only container filesystem.
6. Keep data, credentials, runtime configuration and site-specific launch policy outside the image.
7. Inspect the executable architecture and shared-library dependencies locally.
8. Transfer or publish the Docker/OCI image and convert it to a named SIF image.
9. Inspect the actual library resolution under the appropriate Pawsey Singularity module.
10. Run a representative Slurm job and verify its output before trusting performance results.
11. Review local Docker images, containers and build cache before removing them.

These practices produce an image that is smaller, safer and easier to reuse. More importantly, they separate the stable application environment from the data, configuration and launch policy that vary between users, systems and executions.
