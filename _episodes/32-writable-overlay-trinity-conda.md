---
title: "\"Writable\" containers: bind mounts, overlays, Trinity and Conda"
teaching: 25
exercises: 15

questions:
- Where can a containerised application store files when its container image is read-only?
- How can I make a specific host directory available inside a container?
- Can files created within the container filesystem persist between container runs?
- Can I add software to a container without rebuilding its image?
- How can I avoid storing thousands of individual application or installation files directly on a parallel filesystem?
- What should I use when changes to the container filesystem are needed only during one container run?

objectives:
- Bind mount a host directory at a specific path inside a container
- Create and mount a persistent overlay with a container image
- Use a persistent overlay to store application output and additional software separately from the container image
- Use an overlay to reduce the number of individual files presented to the host filesystem
- Use a temporary writable layer when changes do not need to persist

keypoints:
- Bind mounts make host files directly accessible at selected paths inside a container
- An immutable container image can be combined with a separate writable overlay
- Changes stored in a persistent overlay remain available across container runs
- From the host filesystem's perspective, an overlay consolidates many internal files into a single overlay file
- "`--writable-tmpfs` provides a new disposable writable layer for each container run"
---

### Request an interactive allocation

If you're running this tutorial on a shared system (*e.g.* Setonix at Pawsey), you should use one of the compute nodes rather than the login node. You can do this by requesting an interactive allocation from the scheduler, for instance on Setonix with Slurm (do this if you are not in an `salloc` interactive session yet):

```
$ salloc -N 1 -n 1 -c 4 --reservation=ContainersTraining -t 4:00:00
```
{: .source}

```text
salloc: Granted job allocation 3453895
salloc: Waiting for resource configuration
salloc: Nodes nid000152 are ready for job
```
{: .output}

### Get ready for the hands-on

Before we start, let us ensure we have the files and container image required for the hands-on exercises.

If you have not done so already, move to a suitable working directory and download the following GitHub repository. On Pawsey systems, use your scratch directory; on other HPC or cloud systems, use the equivalent working directory recommended by the system administrators.

```bash
$ cd "$MYSCRATCH"    # On Pawsey systems
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
```
{: .source}

Now move to the working directory for the Trinity example:

```bash
$ cd "$TUTO/demos/trinity"
$ pwd
```
{: .source}

The working directory should be something like:

```text
/path/to/scratch/singularity-containers/demos/trinity
```
{: .output}

Load the singularity module:
```bash
$ module load singularity/4.1.0-nompi
```
{: .source}

Download the Ubuntu container image used throughout this episode and save it in your local image library:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
$ singularity pull "${MY_LOCAL_LIBRARY}/ubuntu--24.04.sif" docker://docker.io/ubuntu:24.04
$ UBUNTU_IMAGE="${MY_LOCAL_LIBRARY}/ubuntu--24.04.sif"
```
{: .source}

### Make a container path writable with a bind mount

As discussed in the basic Singularity episode, bind mounts make host directories available at selected paths inside a container. They are also useful when an application needs to write to a path that belongs to the read-only container filesystem.

A bind mount does not make the directory packaged in the container image writable. Instead, it temporarily covers that path with a host directory. If the host directory is writable, the application can write through the container path, and the resulting files are stored directly in the host directory.

For example, create a host directory and bind mount it at `/run` inside the container:

```bash
$ mkdir -p "$PWD/my_run"
$ singularity exec --bind "$PWD/my_run:/run" "$UBUNTU_IMAGE" \
    bash -c 'touch /run/running-file; ls -lth /run'
```
{: .source}

```text
-rw-r--r-- 1 courses01 courses 0 Sep 17 09:47 running-file
```
{: .output}

The command creates `running-file` in the host directory that was made visible as `/run` inside the container. Back on the host, the file is directly available:

```bash
$ ls -lth my_run
```
{: .source}

```text
-rw-r--r-- 1 courses01 courses 0 Sep 17 09:47 running-file
```
{: .output}

Without the bind mount, `/run` shows the content provided by the container image instead:

```bash
$ singularity exec "$UBUNTU_IMAGE" ls -lth /run
```
{: .source}

```text
drwxr-xr-x 2 root root 32 Sep 11 10:10 systemd
drwxrwxrwt 2 root root  3 Sep 11 10:03 lock
```
{: .output}

The original content of `/run` was not removed or modified by the bind mount. It was hidden for the duration of the previous container run and became visible again when the container was run without that bind mount.


> ## Existing content at the destination path
>
> A bind mount does not merge the host directory with the directory packaged in the container image. For the duration of the container run, the bind-mounted host directory replaces the original directory in the container's visible filesystem. Any content originally present at that path is hidden, but remains unchanged in the container image and becomes visible again when the container is run without that bind mount.
{: .callout}

Bind mounts can also expose host-provided software and libraries inside a container. Later in this training, the MPI episode shows how Pawsey's MPI-enabled Singularity module uses bind mounts, together with library-path and preload settings, to make the host MPI and interconnect environment available to containerised applications.

Bind mounts are generally the simplest approach when files should remain directly accessible on the host. However, every file remains a separate entry on the host filesystem. A persistent overlay can be more suitable when a workflow creates a very large number of small files, or when persistent changes are needed across multiple paths in the container filesystem.

### Create a persistent overlay filesystem

In a previous episode, we saw that the filesystem contained in a SIF container image is **read-only**. In other words, files within that filesystem cannot normally be created, modified or removed.

However, there are situations where, rather than reading and writing files in a bind-mounted host directory, it is useful to store changes persistently within the container filesystem.

A practical use case arises when an application creates a very large number of small files on a host parallel filesystem such as *Lustre*. Storing every file directly on Lustre can consume file quotas and generate substantial workloads for its metadata services.

A persistent overlay can reduce this metadata workload by storing those files within a single overlay file from the host filesystem's perspective. The host filesystem still handles reads and writes to the overlay file, but it does not manage each file inside the overlay as a separate entry.

Singularity supports persistent overlays for this purpose. A persistent overlay stores changes separately from the immutable container image and makes those changes available again whenever the overlay is mounted.

Now let's create a 200 MiB overlay:

```bash
$ singularity overlay create --size 200 my_overlay.ext3
```
{: .source}


### Mount a persistent overlay filesystem

Let's give it a go with the persistent overlay we have just created. We can mount it at container runtime by using the `--overlay` option followed by the overlay filename. By default, the overlay is mounted read-write. We make this explicit by adding the `:rw` suffix to the filename:

```bash
$ singularity shell --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE"
```
{: .source}

Now, every new directory and file that we create from inside the container will be stored in the persistent overlay filesystem.  For instance, from the interactive shell we opened let us try:

```
Singularity> mkdir /australia
Singularity> cd /australia
Singularity> echo perth > wa
Singularity> echo canberra > act
Singularity> exit
```
{: .bash}


> ## Access a pre-existing overlay filesystem
>
> Once exited the container, the newly created directory is not available in the host filesystem.  Try and inspect the content of `/australia` from the host.
>
> > ## Solution (with error)
> >
> > ```
> > $ ls /australia
> > ```
> > {: .bash}
> >
> > ```
> > ls: /australia: No such file or directory
> > ```
> > {: .output}
> {: .solution}
>
> Now try and look for `/australia` from inside a Ubuntu container, *without* mounting the overlay.
>
> > ## Solution (with error)
> >
> > ```
> > $ singularity exec "$UBUNTU_IMAGE" ls /australia
> > ```
> > {: .bash}
> >
> > ```
> > ls: /australia: No such file or directory
> > ```
> > {: .output}
> {: .solution}
>
> But, if we run another container and mount the overlay, the files will still be there.  Try and `ls` from inside a container, using the appropriate flag.
>
> > ## Solution
> >
> > ```
> > $ singularity exec --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE" ls /australia
> > ```
> > {: .bash}
> >
> > ```
> > act wa
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


The newly created directories and files persist in the overlay and can be accessed again in future container runs whenever `my_overlay.ext3` is mounted. Data files may also be accessed with other compatible container images, although software stored in an overlay may depend on the original container image.


### Run a Trinity genome assembly from inside the container

A subdirectory in the directory we are in, `trinity_test_data/`, contains sample inputs for a genome assembly, coming from the `Docker/test_data/` subset in the [Trinity Github repo](https://github.com/trinityrnaseq/trinityrnaseq).


> ## Create the output directory within an OverlayFS
>
> For this exercise we are going to reuse the overlay file `my_overlay.ext3`.
> To begin with, use the Ubuntu container image `ubuntu--24.04.sif` to create the directory `/trinity_out_dir` in the OverlayFS.
>
> > ## Solution
> >
> > ```
> > $ singularity exec --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE" mkdir /trinity_out_dir
> > $ singularity exec --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE" ls -ltd /trinity_out_dir
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


Now, let's download the Trinity container image from Docker Hub, `trinityrnaseq/trinityrnaseq:2.8.6`:

```
$ singularity pull "${MY_LOCAL_LIBRARY}/trinityrnaseq--2.8.6.sif" docker://docker.io/trinityrnaseq/trinityrnaseq:2.8.6
```
{: .bash}

Now, we're going to run a test assembly with our sample dataset, using the directory we just created in the OverlayFS to write the outputs.  In this small case, only about a hundred output files are created.  However, the Trinity assembly workflow can easily produce up to a million files, making the use of OverlayFS an interesting choice to reduce the workload on a host parallel filesystem.


> ## Run the Trinity workflow
>
> This is the command to be run:
>
> ```
> $ Trinity \
>     --seqType fq --left trinity_test_data/reads.left.fq.gz \
>     --right trinity_test_data/reads.right.fq.gz \
>     --max_memory 1G --CPU 1 --output <OUTPUT-DIRECTORY>
> ```
> {: .bash}
>
> Can you run it using Singularity and OverlayFS?  (**Hint**: you will need to specify the appropriate directory for `--output`)
> Expect the run to take approximately 2-3 minutes.
>
> > ## Solution
> >
> > ```
> > $ TRINITY_IMAGE="${MY_LOCAL_LIBRARY}/trinityrnaseq--2.8.6.sif"
> > $ ls -lth "$TRINITY_IMAGE"
> > $ singularity exec --overlay "my_overlay.ext3:rw" "$TRINITY_IMAGE" \
> >     Trinity \
> >     --seqType fq --left trinity_test_data/reads.left.fq.gz \
> >     --right trinity_test_data/reads.right.fq.gz \
> >     --max_memory 1G --CPU 1 --output /trinity_out_dir
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


All of our outputs are stored in the OverlayFS, so we need to use a Singularity container to inspect them:

```bash
$ singularity exec --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE" ls /trinity_out_dir
```
{: .source}

```text
Trinity.fasta		      inchworm.K25.L25.DS.fa.finished  pipeliner.1217836.cmds
Trinity.fasta.gene_trans_map  inchworm.kmer_count	       read_partitions
Trinity.timing		      insilico_read_normalization      recursive_trinity.cmds
both.fa			      jellyfish.kmers.fa	       recursive_trinity.cmds.completed
both.fa.ok		      jellyfish.kmers.fa.histo	       recursive_trinity.cmds.ok
both.fa.read_count	      left.fa.ok		       right.fa.ok
chrysalis		      partitioned_reads.files.list     scaffolding_entries.sam
inchworm.K25.L25.DS.fa	      partitioned_reads.files.list.ok
```
{: .output}

Now let's copy the assembled sequence and transcripts, `Trinity.fasta*`, in the current directory:

```bash
$ singularity exec --overlay "my_overlay.ext3:rw" "$UBUNTU_IMAGE" bash -c 'cp -p /trinity_out_dir/Trinity.fasta* ./'
```
{: .source}

Note how we're wrapping the copy command within `bash -c`; this is to defer the evaluation of the `*` wildcard to when the container runs the command.

We've run the entire workflow within the OverlayFS, and got only the two relevant output files out in the host filesystem!

```
$ ls -l Trinity.fasta*
```
{: .bash}

```
-rw-r--r-- 1 courses01 courses 171507 Nov  4 05:49 Trinity.fasta
-rw-r--r-- 1 courses01 courses   2818 Nov  4 05:49 Trinity.fasta.gene_trans_map
```
{: .output}


### Installing software with Conda/Mamba inside an overlay

Persistent overlays can be used to extend an existing container without modifying the container image itself. The container image remains immutable, while additional software is installed separately in the overlay and becomes available whenever that overlay is mounted. In practice, you would normally choose a container image that already provides most of the software stack you need, then use an overlay to add any missing packages or tools.

This approach is particularly useful for Conda/Mamba installations, which typically consist of many thousands of small files. Storing those files directly on a host parallel filesystem can consume file quotas and create substantial metadata workloads. Installing the Conda/Mamba environment inside an overlay consolidates all those installation files into one overlay file from the host filesystem's perspective, while still making the installed software available inside the container.

Let's work in another directory:

```
$ cd "$TUTO/demos/conda_overlay"
```
{: .bash}

And let's create a new overlay for this example. This time, we'll make it considerably larger, since a Conda/Mamba installation can easily take up a few gigabytes (note that the `--size` value in the following command is specified in MiB):

```
$ singularity overlay create --size 5000 my_conda_overlay.ext3
```
{: .bash}

> ## Choose the overlay size carefully
>
> Singularity overlays have a fixed size that must be chosen at creation time. If you do not reserve enough space, the installation will fail partway through and you will need to create a larger overlay and start again. When in doubt, err on the generous side.
>
> An alternative approach for storing large collections of small files is to package them into a SquashFS file (see the Pawsey's documentation page: ["How to use SquashFS to avoid file quota issues"](https://pawsey.atlassian.net/wiki/spaces/US/pages/51927678/How+to+use+SquashFS+to+avoid+file+quota+issues)). Unlike overlays, SquashFS files do not require preallocating storage space. As SquashFS is a filesystem packaging technology rather than a container technology, it is outside the scope of this lesson.
{: .callout}

We're going to use `ubuntu--24.04.sif` again for this example. This minimal Ubuntu container image does not ship with `wget` or `curl`, so rather than downloading the Miniforge installer *from inside* the container, we'll download it first on the **host** into the current directory. Singularity will bind mount the current working directory into the container by default, making the installer available from inside the container later. On the x86-64 system used in this lesson, download the corresponding Miniforge installer:

```bash
$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
```
{: .source}

Now, let's open a shell in the container with our new overlay mounted read-write:

```bash
$ singularity shell --overlay "my_conda_overlay.ext3:rw" "$UBUNTU_IMAGE"
```
{: .source}

From inside the container, run the Miniforge installer and use `-p /opt/conda` to select the installation directory. Because `/opt/conda` is part of the container filesystem and does not exist in the read-only container image, the new directory and its contents will be stored in the overlay.

Do **not** create `/opt/conda` beforehand, as the installer expects to create the installation directory itself.

```bash
Singularity> bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda
```
{: .source}
 
The `-b` flag runs the installer in batch, or non-interactive, mode. The `-p /opt/conda` option installs Miniforge under `/opt/conda` instead of using a location under `$HOME`. The installation is stored in the overlay rather than being added to the immutable container image. From the host filesystem's perspective, its files remain contained within the overlay file.

Once Miniforge is installed, add its `bin` directory to `PATH` so that its commands are available in the current shell:

```bash
Singularity> export PATH="/opt/conda/bin:$PATH"
```
{: .source}

This change applies only to the current shell. It will be lost when we exit, so `/opt/conda/bin` must be added to `PATH` again in later container runs that need to use the installed commands.

From this point, you can use `conda` and `mamba` as usual, including creating and managing environments. If you need to activate a Conda environment in the current shell, first source `/opt/conda/etc/profile.d/conda.sh`. For this example, no environment activation is needed.

Let's use `mamba` to install the *BWA* aligner from the *bioconda* channel:

```bash
Singularity> mamba install -y -c conda-forge -c bioconda bwa
Singularity> mamba clean --all --yes
Singularity> exit
```
{: .source}

The `mamba clean` command removes package caches that are no longer needed, freeing some space in the overlay. After exiting the container, we can remove the Miniforge installer from the host:

```bash
$ rm Miniforge3-Linux-x86_64.sh
```
{: .source}

> ## Use the installed software in another container run
>
> Exit the interactive shell, then use `singularity exec` with the same container image and overlay to check that `bwa` remains available.
>
> Remember that the `PATH` change made in the previous interactive shell does not persist after exiting. Add `/opt/conda/bin` to `PATH` inside this new container run before executing `bwa`.
>
> > ## Solution
> >
> > ```bash
> > $ singularity exec --overlay "my_conda_overlay.ext3:rw" "$UBUNTU_IMAGE" \
> >     bash -c 'export PATH=/opt/conda/bin:$PATH && bwa'
> > ```
> > {: .source}
> >
> > ```text
> > Program: bwa (alignment via Burrows-Wheeler transformation)
> > Version: 0.7.17-r1188
> > Contact: Heng Li <lh3@sanger.harvard.edu>
> > ...
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

The container image has not been modified. The BWA installation persists because it is stored in `my_conda_overlay.ext3` and becomes available whenever that overlay is mounted. However, shell settings such as the modified `PATH` do not persist between container runs, so they must be set again in each new shell or command.

The same overlay may also be mounted with other compatible container images. However, software installed in the overlay may depend on libraries or other components provided by the original container image. For reliable and reproducible use, mount the overlay with the same container image unless you have verified compatibility with another one.

### Ephemeral writable containers

Sometimes an application needs to create temporary files within the container filesystem, but you do not need those changes to persist after the run. Examples include software that creates cache files, lock files, runtime state under `/run`, temporary data under `/tmp`, or configuration files that are only needed while the application is running.

In these situations, creating a persistent overlay may be unnecessary. Instead, Singularity can make the normally read-only parts of a SIF container filesystem temporarily writable using the `--writable-tmpfs` flag. Those changes are stored in a temporary writable layer and are discarded when the container run ends. Writes to bind-mounted host directories are not stored in this temporary layer.

For example, let's start an interactive shell with a temporary writable layer:

```bash
$ singularity shell --writable-tmpfs "$UBUNTU_IMAGE"
```
{: .source}

Inside the container, create a file at the root of the container filesystem and confirm that it exists:

```bash
Singularity> touch /temporary-file
Singularity> ls -l /temporary-file
```
{: .source}

```text
-rw-r--r-- 1 user group 0 Sep 17 09:00 /temporary-file
```
{: .output}

Exit the container:

```bash
Singularity> exit
```
{: .source}

The temporary writable layer is discarded when the container run ends. To confirm this, start another container run with a new temporary writable layer and check for the file:

```bash
$ singularity exec --writable-tmpfs "$UBUNTU_IMAGE" ls /temporary-file
```
{: .source}

```text
ls: cannot access '/temporary-file': No such file or directory
```
{: .output}

The file existed during the previous interactive shell, but it was discarded along with that shell's temporary writable layer. Each use of `--writable-tmpfs` creates a new temporary layer; changes from an earlier container run are not recovered in subsequent runs.

Use `--writable-tmpfs` for temporary, disposable changes. Use a persistent overlay when changes need to survive across container runs, and use a bind mount when files should remain directly accessible on the host filesystem.
