---
title: "Set up writable containers: another bio example with Trinity"
teaching: 10
exercises: 10
questions:
objectives:
- Create and mount a writable overlay filesystem
- Make a container temporarily writable
keypoints:
- Use linux tools from a Ubuntu container to create a filesystem image file
- Mount a filesystem image using the flag `--overlay`
- Make a SIF container ephemerally writable with the flag `--writable-tmpfs`
---

### Request an interactive allocation

If you're running this tutorial on a shared system (*e.g.* Setonix at Pawsey), you should use one of the compute nodes rather than the login node. You can do this by requesting an interactive allocation from the scheduler, for instance on Setonix with Slurm (do this if you are not in an `salloc` interactive session yet):

```
$ salloc -N 1 -n 1 -c 16 --reservation=ContainersTraining -t 4:00:00
```
{: .source}

```text
salloc: Granted job allocation 3453895
salloc: Waiting for resource configuration
salloc: Nodes nid000152 are ready for job
```
{: .output}

### Get ready for the hands-on

Before we start, let us ensure we have got the required files to run the tutorials.

If you haven't done so already, move to a suitable working directory and download the following GitHub repository. On Pawsey systems, use your scratch directory; on other HPC or cloud systems, use the equivalent working directory recommended by the system administrators.

```bash
$ cd "$MYSCRATCH"    # On Pawsey systems
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
```
{: .source}

Now `cd` to the working directory. In this case:
```bash
$ cd $TUTO/demos/trinity
$ pwd
```

The working directory should be something like:
```text
/path/to/scratch/singularity-containers/demos/trinity
```
{: .source}


### Create a persistent overlay filesystem

In a previous episode, we've seen that Singularity containers are **read-only**.  In other words, you cannot create or edit any files inside their filesystem.

However, there can be instances where, rather than reading/writing files in the host filesystem, it would instead come handy to persistently store them inside the container filesystem.
A practical user case is when using a host parallel filesystem such as *Lustre* to run applications that create a large number (*e.g.* millions) of small files.  This practice creates a huge workload on the metadata servers of the filesystem, degrading its performance.  In this context, significant performance benefits can be achieved by reading/writing these files inside the container.

Singularity offers a feature to achieve this, called *OverlayFS* and has a dedicated syntax that can be used together with ubuntu images to create these files.

Now let's create the overlay file with `SIZE=200MB`:

```bash
$ module load singularity/4.1.0-nompi
$ export SIZE="200"
$ export FILE="my_overlay.ext3"
$ singularity overlay create --size $SIZE $FILE
```
{: .source}


### Mount a persistent overlay filesystem

Let's give it a go with the persistent filesystem image we have just created.  We can mount it at container runtime by using the flag `--overlay` followed by the image filename:

First of all, let's download the ubuntu image to use and save it in our own image library:

```bash
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
$ singularity pull "${MY_LOCAL_LIBRARY}/ubuntu--24.04.sif" docker://docker.io/ubuntu:24.04
```
{: .source}

```bash
$ UBUNTU_IMAGE="${MY_LOCAL_LIBRARY}/ubuntu--24.04.sif"
$ singularity shell --overlay my_overlay.ext3 $UBUNTU_IMAGE
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
> > ## Solution
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
> > ## Solution
> >
> > ```
> > $ singularity exec $SINGULARITY_IMAGE ls /australia
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
> > $ singularity exec --overlay my_overlay $SINGULARITY_IMAGE ls /australia
> > ```
> > {: .bash}
> >
> > ```
> > act wa
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


Note how the newly created directories and files are persistent, therefore can be re-accessed and re-used in future runs, even by containers instantiated from different images.  All we have to do is to mount the filesystem image `my_overlay.ext3`.


### Run a Trinity genome assembly from inside the container

A subdirectory in the directory we are in, `trinity_test_data/`, contains sample inputs for a genome assembly, coming from the `Docker/test_data/` subset in the [Trinity Github repo](https://github.com/trinityrnaseq/trinityrnaseq).


> ## Create the output directory within an OverlayFS
>
> For this exercise we are going to reuse the filesystem image file `my_overlay.ext3`.
> To begin with, use the Singularity image `ubuntu--24.04.sif` to create the directory `/trinity_out_dir` in the OverlayFS.
>
> > ## Solution
> >
> > ```
> > $ singularity exec --overlay my_overlay.ext3 $UBUNTU_IMAGE mkdir /trinity_out_dir
> > $ singularity exec --overlay my_overlay.ext3 $UBUNTU_IMAGE ls -ltd /trinity_out_dir
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


Now, let's download the Trinity image from Docker hub, `trinityrnaseq/trinityrnaseq:2.8.6`:

```
$ singularity pull docker://trinityrnaseq/trinityrnaseq:2.15.2
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
> > $ singularity exec --overlay my_overlay.ext3 "$TRINITY_IMAGE" \
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
$ singularity exec --overlay my_overlay.ext3 "$UBUNTU_IMAGE" ls /trinity_out_dir
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
$ singularity exec --overlay my_overlay.ext3 "$UBUNTU_IMAGE" bash -c 'cp -p /trinity_out_dir/Trinity.fasta* ./'
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

Writing large numbers of output files isn't the only scenario where a persistent overlay is useful.  Another common case is *installing software*: package managers such as *Conda* (and its faster drop-in replacement, *Mamba*) typically create installations made up of many thousands of small files.  Installing directly on a host parallel filesystem can therefore run into the same file quota and metadata performance problems we discussed above.  We can instead install the whole Conda/Mamba environment *inside* a persistent overlay, keeping all of those small files neatly packed away in a single image file on the host filesystem.

Let's create a new, empty directory to work in, and a fresh overlay image dedicated to this example.  This time, we'll make it considerably bigger than `my_overlay.ext3`, since a Conda/Mamba installation can easily take up a few gigabytes (note that size is in MB):

```
$ cd $TUTO/demos/conda_overlay
$ export SIZE="5000"
$ export FILE="my_conda_overlay.ext3"
$ singularity overlay create --size $SIZE $FILE
$ ls -lat
```
{: .bash}

> ## Mind the size!
>
> Singularity overlays have a fixed size that must be chosen at creation time. If you do not reserve enough space, the installation will fail partway through and you will need to create a larger overlay and start again. When in doubt, err on the generous side.
>
> An alternative approach for storing large collections of small files is to package them into a SquashFS file (see the [Pawsey documentation on SquashFS](https://pawsey.atlassian.net/wiki/spaces/US+to+avoid+file+quota+issues). Unlike overlays, SquashFS files do not require preallocating storage space. As SquashFS is a filesystem packaging technology rather than a container technology, it is outside the scope of this lesson.
{: .callout}

We're going to use again `ubuntu--24.04.sif` image for this.  Such a small image doesn't ship with `wget` or `curl`, so rather than downloading the miniconda installer *from inside* the container, we download it first directly with **host** command into our current directory (which Singularity will bind mount later by default when using the container):

```bash
$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
```
{: .source}

> ## Which base image?
>
> Any Linux container image will do here: the [Miniforge](https://github.com/conda-forge/miniforge) installer is self-contained and only needs `bash` and core utilities, and Conda/Mamba then bring their own Python and SSL libraries.  We stick with `ubuntu--24.04.sif` for consistency with the rest of this episode;
{: .callout}

Now, let's open a shell into the container, mounting our new overlay as read-write with the `:rw` suffix (this is the default, but it doesn't hurt to be explicit):

```bash
$ singularity shell --overlay $FILE:rw "$UBUNTU_IMAGE"
```
{: .source}

From inside the container, we run the installer, pointing it at a path at the root of the overlay with `-p`.  Note we do **not** create `/opt/conda` beforehand: the installer refuses to write into a directory that already exists.

```
Singularity> bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda
```
{: .bash}

The `-b` flag runs the installer in batch (non-interactive) mode, and `-p /opt/conda` installs into the overlay rather than under `$HOME`.  All the files land inside the overlay image, not on the host filesystem.

Once installed, we can activate the environment and use `mamba` to install whatever packages we need.  Let's add the *bwa* aligner, which lives on the *bioconda* channel:

```
Singularity> source /opt/conda/etc/profile.d/conda.sh
Singularity> export PATH=/opt/conda/bin:$PATH
Singularity> mamba install -y -c bioconda -c conda-forge bwa
Singularity> mamba clean --all --yes
Singularity> exit
```
{: .bash}

We ran `mamba clean` before exiting to remove downloaded package caches we no longer need, saving some space in the overlay.  Back on the host, we can delete the installer script:

```
$ rm Miniforge3-Linux-x86_64.sh
```
{: .bash}

> ## Use the installed software from a fresh container
>
> Once exited, start a **new** `singularity exec` from the *same* Ubuntu image, mounting `my_conda_overlay.ext3` again, and check that `bwa` is available and runs.  (**Hint**: you will need to add `/opt/conda/bin` to the `PATH` *inside* the container before calling `bwa`).
>
> > ## Solution
> >
> > ```bash
> > $ singularity exec --overlay my_conda_overlay.ext3 "$UBUNTU_IMAGE" bash -c 'export PATH=/opt/conda/bin:$PATH && bwa'
> > ```
> > {: .source}
> >
> > ```
> > Program: bwa (alignment via Burrows-Wheeler transformation)
> > Version: 0.7.17-r1188
> > Contact: Heng Li <lh3@sanger.harvard.edu>
> > ...
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Just like with the Trinity example earlier, the software we installed persists inside the overlay image file, and can be reused across different container runs, even from different container images, simply by mounting `my_conda_overlay.ext3` again with `--overlay`.  This makes overlays a handy way to keep bulky, many-file software installations off your quota-limited host filesystem, while still being able to bring them along with any Singularity container you like.


### Ephemeral writable containers

Sometimes an application needs to create temporary files within the container filesystem, but you do not need those changes to persist after the run. Examples include software that creates cache files, lock files, runtime state under `/run`, temporary data under `/tmp`, or configuration files that are only needed while the application is running.

In these situations, creating a persistent overlay may be unnecessary. Instead, Singularity can make a SIF container temporarily writable using the `--writable-tmpfs` flag. Any changes made to the container filesystem are stored in a temporary writable layer and are discarded when the container exits.

For example, let's create a file under `/tmp`:

```bash
$ singularity exec --writable-tmpfs "$UBUNTU_IMAGE" touch /tmp/temporary-file
```
{: .source}

The command succeeds because the container filesystem is writable for the duration of the run.

To verify that the change did not persist, launch a fresh container and check for the file:

```bash
$ singularity exec "$UBUNTU_IMAGE" ls /tmp/temporary-file
```
{: .source}

The output should be:

```text
ls: cannot access '/tmp/temporary-file': No such file or directory
```
{: .output}

The file disappeared because it was created in the temporary writable layer provided by `--writable-tmpfs`.

#### When `--writable-tmpfs` is not enough

Some applications need more writable space than is available in the temporary writable layer, or need their changes to persist across multiple container runs. In these situations, use a persistent overlay instead.

Another common approach is to bind mount a host directory at the location where the application needs write access. For example:

```bash
$ mkdir -p "$PWD/my_run"
$ singularity exec --bind "$PWD/my_run:/run" "$UBUNTU_IMAGE" touch /run/running-file
```
{: .source}

The file is written into the host directory `my_run` and therefore remains available after the container exits:

```bash
$ ls my_run
```
{: .source}

```text
running-file
```
{: .output}

Use `--writable-tmpfs` for temporary, disposable changes, and use persistent overlays when changes need to survive across container runs.
