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


### Create a persistent overlay filesystem

In a previous episode, we've seen that Singularity containers are **read-only**.  In other words, you cannot create or edit any files inside their filesystem.

However, there can be instances where, rather than reading/writing files in the host filesystem, it would instead come handy to persistently store them inside the container filesystem.  
A practical user case is when using a host parallel filesystem such as *Lustre* to run applications that create a large number (*e.g.* millions) of small files.  This practice creates a huge workload on the metadata servers of the filesystem, degrading its performance.  In this context, significant performance benefits can be achieved by reading/writing these files inside the container.

Singularity offers a feature to achieve this, called *OverlayFS*.

Let us cd into `demos/trinity`:

```
$ cd $TUTO/demos/trinity
```
{: .bash}

and then discuss how to use the Linux tools `dd` and `mkfs.ext3` to create and format an empty *ext3* filesystem image, which we will call `my_overlay`.  These Linux tools typically require `sudo` privileges to run.  However, we can bypass this requirement by using the ones provided inside a standard *Ubuntu* container.  The following command looks a bit cumbersome, but is indeed just an idiomatic syntax to achieve our goal with Singularity (up to versions 3.7.x):

```
$ export COUNT="200"
$ export BS="1M"
$ export FILE="my_overlay"
$ singularity exec docker://ubuntu:18.04 bash -c " \
    mkdir -p overlay_tmp/upper overlay_tmp/work && \
    dd if=/dev/zero of=$FILE count=$COUNT bs=$BS && \
    mkfs.ext3 -d overlay_tmp $FILE && \
    rm -rf overlay_tmp \
    "
```
{: .bash}

```
200+0 records in
200+0 records out
209715200 bytes (210 MB, 200 MiB) copied, 0.362418 s, 579 MB/s
mke2fs 1.44.1 (24-Mar-2018)
ext2fs_check_if_mount: Can't check if filesystem is mounted due to missing mtab file while determining whether my_overlay is mounted.
Discarding device blocks: done
Creating filesystem with 204800 1k blocks and 51200 inodes
Filesystem UUID: b14f9b2a-188d-4c19-8e6d-38a568f6efe1
Superblock backups stored on blocks:
	8193, 24577, 40961, 57345, 73729

Allocating group tables: done
Writing inode tables: done
Creating journal (4096 blocks): done
Copying files into the device: done
Writing superblocks and filesystem accounting information: done

```
{: .output}

Here we have wrapped four commands into a single bash call from a container, just for the convenience of running it once.  We've also defined shell variables for better clarity.  
What are the single commands doing?  
We are creating (and then deleting at the end) two service directories, `overlay_tmp/upper` and `overlay_tmp/work`, that will be used by the command `mkfs.ext3`.  
The `dd` command creates a file named `my_overlay`, made up of blocks of zeros, namely with `count` blocks of size `bs` (the unit here is *megabytes*); the product `count*bs` gives the total file size in bytes, in this case corresponding to *200 MB*.
The command `mkfs.ext3` is then used to format the file as a *ext3* filesystem image, that will be usable by Singularity.  Here we are using the service directory we created, `my_overlay`, with the flag `-d`, to tell `mkfs` we want the filesystem to be owned by the same owner of this directory, *i.e.* by the current user.  If we skipped this option, we would end up with a filesystem that is writable only by *root*, not very useful.

Note how, starting from version 3.8, Singularity offers a dedicated syntax that wraps arounds the commands above, providing a simpler interface (here size must be in MB):

```
$ export SIZE="200"
$ export FILE="my_overlay"
$ singularity overlay create --size $SIZE $FILE
```
{: .bash}


### Mount a persistent overlay filesystem

Let's give it a go with the persistent filesystem image we have just created.  We can mount it at container runtime by using the flag `--overlay` followed by the image filename:

```
$ singularity shell --overlay my_overlay docker://ubuntu:18.04
```
{: .bash}

Now, every new directory and file that we create from inside the container will be stored in the persistent overlay filesystem.  For instance, from the interactive shell we opened let us try:

```
Singularity> mkdir /australia
Singularity> cd /australia
Singularity> echo perth >wa
Singularity> echo canberra >act
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
> > $ singularity exec docker://ubuntu:18.04 ls /australia
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
> > $ singularity exec --overlay my_overlay docker://ubuntu:18.04 ls /australia
> > ```
> > {: .bash}
> >
> > ```
> > act wa
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


Note how the newly created directories and files are persistent, therefore can be re-accessed and re-used in future runs, even by containers instantiated from different images.  All we have to do is to mount the filesystem image `my_overlay`.


### Run a Trinity genome assembly from inside the container

A subdirectory in the directory we are in, `trinity_test_data/`, contains sample inputs for a genome assembly, coming from the `Docker/test_data/` subset in the [Trinity Github repo](https://github.com/trinityrnaseq/trinityrnaseq).


> ## Create the output directory within an OverlayFS
>
> For this exercise we are going to reuse the filesystem image file `my_overlay`.  
> To begin with, use the Singularity image `ubuntu:18.04` to create the directory `/trinity_out_dir` in the OverlayFS.
>
> > ## Solution
> >
> > ```
> > $ singularity exec --overlay my_overlay docker://ubuntu:18.04 mkdir /trinity_out_dir
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


Now, let's download the Trinity image from Docker hub, `trinityrnaseq/trinityrnaseq:2.8.6`:

```
$ singularity pull docker://trinityrnaseq/trinityrnaseq:2.8.6
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
> > $ singularity exec --overlay my_overlay trinityrnaseq_2.8.6.sif \
> >     Trinity \
> >     --seqType fq --left trinity_test_data/reads.left.fq.gz \
> >     --right trinity_test_data/reads.right.fq.gz \
> >     --max_memory 1G --CPU 1 --output /trinity_out_dir
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


All of our outputs are stored in the OverlayFS, so we need to use a Singularity container to inspect them:

```
$ singularity exec --overlay my_overlay docker://ubuntu:18.04 ls /trinity_out_dir
```
{: .bash}

```
Trinity.fasta		      both.fa.read_count	       insilico_read_normalization   partitioned_reads.files.list.ok   recursive_trinity.cmds.ok
Trinity.fasta.gene_trans_map  chrysalis			       jellyfish.kmers.fa	     pipeliner.18881.cmds	       right.fa.ok
Trinity.timing		      inchworm.K25.L25.DS.fa	       jellyfish.kmers.fa.histo      read_partitions		       scaffolding_entries.sam
both.fa			      inchworm.K25.L25.DS.fa.finished  left.fa.ok		     recursive_trinity.cmds
both.fa.ok		      inchworm.kmer_count	       partitioned_reads.files.list  recursive_trinity.cmds.completed
```
{: .output}

Now let's copy the assembled sequence and transcripts, `Trinity.fasta*`, in the current directory:

```
$ singularity exec --overlay my_overlay docker://ubuntu:18.04 bash -c 'cp -p /trinity_out_dir/Trinity.fasta* ./'
```
{: .bash}

Note how we're wrapping the copy command within `bash -c`; this is to defer the evaluation of the `*` wildcard to when the container runs the command.

We've run the entire workflow within the OverlayFS, and got only the two relevant output files out in the host filesystem!

```
$ ls -l Trinity.fasta*
```
{: .bash}

```
-rw-r--r-- 1 ubuntu ubuntu 171507 Nov  4 05:49 Trinity.fasta
-rw-r--r-- 1 ubuntu ubuntu   2818 Nov  4 05:49 Trinity.fasta.gene_trans_map
```
{: .output}


### Installing software with Conda/Mamba inside an overlay

Writing large numbers of output files isn't the only scenario where a persistent overlay is useful.  Another common case is *installing software*: package managers such as *Conda* (and its faster drop-in replacement, *Mamba*) typically create installations made up of many thousands of small files.  Installing directly on a host parallel filesystem can therefore run into the same file quota and metadata performance problems we discussed above.  We can instead install the whole Conda/Mamba environment *inside* a persistent overlay, keeping all of those small files neatly packed away in a single image file on the host filesystem.

Let's create a new, empty directory to work in, and a fresh overlay image dedicated to this example.  This time, we'll make it considerably bigger than `my_overlay`, since a Conda/Mamba installation can easily take up a few gigabytes:

```
$ mkdir -p $TUTO/demos/conda_overlay
$ cd $TUTO/demos/conda_overlay
$ export SIZE="5000"
$ export FILE="my_conda_overlay"
$ singularity overlay create --size $SIZE $FILE
```
{: .bash}

> ## Mind the size!
>
> Unlike the SquashFS images we could have used instead (see the [Pawsey documentation on SquashFS](https://pawsey.atlassian.net/wiki/spaces/US/pages/51927678/How+to+use+SquashFS+to+avoid+file+quota+issues) for that alternative approach), a Singularity overlay has its size fixed at creation time.  If you don't reserve enough space, the installation will fail partway through and you will have to create a bigger overlay and start again.  When in doubt, err on the generous side.
{: .callout}

We're going to use a minimal `ubuntu:18.04` container for this.  Such a small base image doesn't ship with `wget` or `curl`, so rather than downloading the installer *from inside* the container, we download it first on the **host**, in our current directory (which Singularity bind mounts into the container by default):

```
$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
```
{: .bash}

> ## Which base image?
>
> Any Linux container image will do here: the [Miniforge](https://github.com/conda-forge/miniforge) installer is self-contained and only needs `bash` and core utilities, and Conda/Mamba then bring their own Python and SSL libraries.  We stick with `ubuntu:18.04` for consistency with the rest of this episode; if you'd rather fetch the installer from within the container, pick an image that already includes `curl` or `wget` (for instance a Pawsey base image such as `docker://quay.io/pawsey/mpich-base`).
{: .callout}

Now, let's open a shell into the container, mounting our new overlay as read-write with the `:rw` suffix (this is the default, but it doesn't hurt to be explicit):

```
$ singularity shell --overlay $FILE:rw docker://ubuntu:18.04
```
{: .bash}

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
> Once exited, start a **new** `singularity exec` from the *same* Ubuntu image, mounting `my_conda_overlay` again, and check that `bwa` is available and runs.  (**Hint**: you will need to add `/opt/conda/bin` to the `PATH` *inside* the container before calling `bwa`).
>
> > ## Solution
> >
> > ```
> > $ singularity exec --overlay my_conda_overlay docker://ubuntu:18.04 bash -c 'export PATH=/opt/conda/bin:$PATH && bwa'
> > ```
> > {: .bash}
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

Just like with the Trinity example earlier, the software we installed persists inside the overlay image file, and can be reused across different container runs, even from different container images, simply by mounting `my_conda_overlay` again with `--overlay`.  This makes overlays a handy way to keep bulky, many-file software installations off your quota-limited host filesystem, while still being able to bring them along with any Singularity container you like.


### Ephemeral writable containers

In some situations, you might need your container to be writable not to store persistent output files, but just to write temporary service files.  
*E.g.* this can happen with applications that want to write a dot-file in your home, such as a Python package, or containerised Jupyter notebooks that need to write runtime information under `/run`.  
In this context, a persistent overlay filesystem might require more work than is desired.  There are alternative, simpler ways to set this up.

Singularity has a flag for rendering containers from SIF image files ephemerally writable.  `--writable-tmpfs` will allocate a small amount of RAM for this purpose (configured by the sys admins, by default just a bunch of MB), *e.g.*:

```
$ singularity exec --writable-tmpfs docker://ubuntu:18.04 touch ~/write-to-home
```
{: .bash}

Unless `$HOME` is bind mounted to the container (for security reasons it shouldn't be), the newly created file will be gone after the container exits:

```
$ ls ~/write-to-home
```
{: .bash}

```
ls: /home/ubuntu/write-to-home: No such file or directory
```
{: .output}

There are situations where `--writable-tmpfs` is not usable, in particular if you are trying to write to a directory owned by *root*, such as `/run`.  
In this case, the solution is to create a host directory and bind mount it as the path you need to write into, *e.g.*:

```
$ mkdir ~/my_run
$ SINGULARITY_BINDPATH="~/my_run:/run,$SINGULARITY_BINDPATH"
$ singularity exec docker://ubuntu:18.04 touch /run/running-file
```
{: .bash}

In this case, the file will also persist in the host directory after the container exits.
