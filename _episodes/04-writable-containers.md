---
title: "Set up writable containers: another bioinformatics example with Trinity"
teaching: 10
exercises: 10
questions:
objectives:
- Learn how to create and mount a writable overlay filesystem
- Learn how to make a container temporarily writable
keypoints:
- "Use linux tool from a Ubuntu container to create a filesystem image file"
- "Mount a filesystem image using the flag `--overlay`"
- "Make a SIF container ephemerally writable with the flag `--writable-tmpfs`"
---


### Create a persistent overlay filesystem

There can be instances where, rather than reading/writing files in the host filesystem, it would instead come handy to be able to persistently store them inside the container filesystem.  
A practical user case is when using a host parallel filesystem such as *Lustre* to run applications that create a large number (*e.g.* millions) of small files. This practice creates a huge workload on the metadata servers of the filesystem, degrading its performance. In this context, significant performance benefits can be achieved by reading/writing these files inside the container.

Singularity offers a feature to achieve this, called *OverlayFS*.

Let us cd into `demos/04_trinity`:

```
$ cd $SC19/demos/04_trinity
```
{: .bash}

and then discuss how to use the Linux tools `dd` and `mkfs.ext3` to create and format an empty ext3 file system image, which we will call `my_overlay`. These tools typically require `sudo` privileges to run. However, we can bypass this requirement by using the ones provided inside a standard *Ubuntu* container. The following command looks a bit cumbersome, but is indeed just an idiomatic syntax to achieve our goal:

```
singularity exec library://ubuntu:18.04 bash -c " \
  mkdir -p overlay_tmp/upper && \
  dd if=/dev/zero of=my_overlay count=1048576 bs=1024 && \
  mkfs.ext3 -d overlay_tmp my_overlay && \
  rm -rf overlay_tmp \
  "
```
{: .bash}

```

```
{: .output}

Here we have wrapped four commands into a single bash call from a container, just for the convenience of running it once. What are the single commands doing?  
We are creating (and then deleting at the end) a service directory, `overlay_tmp/upper`, that will be used by the command `mkfs.ext3`.  
The `dd` command creates a file named `my_overlay`, made up of blocks of zeros, namely with `count` blocks of size `bs` (the unit here is *bytes*); the product `count*bs` gives the total file size in bytes, in this case corresponding to *1 GB*.  
The command `mkfs.ext3` is then used to format the file as a *ext3* filesystem image, that will be usable by Singularity. Here we are using the service directory we created, `my_overlay`, with the flag `-d`, to tell `mkfs` we want the filesystem to be owned by the same owner of this directory, *i.e.* by the current user. If we skipped this option, we would end up with a filesystem that is writable only by *root*, not very useful.


### Mount a persistent overlay filesystem

Let's give it a go with the persistent filesystem image we have just created. We can mount it at container runtime by using the flag `--overlay` followed by the image filename: 

```
$ singularity shell --overlay my_overlay library://ubuntu:18.04
```
{: .bash}

Now, every new directory and file that we create from inside the container will be stored in the persistent overlay filesystem. For instance, from the interactive shell we opened let us try:

```
$ mkdir /australia
$ cd /australia
$ echo perth >wa
$ echo canberra >act
$ exit
```
{: .bash}


> ## Access a pre-existing overlay filesystem
> 
> Once exited the container, the newly created directory is not available in the host filesystem. Try and inspect the content of `/australia` from the host.
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
> But, if we run another container and mount the overlay, the files will still be there. Try and `ls` from inside a container, using the appropriate flag.
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec --overlay my_overlay library://ubuntu:18.04 ls /australia
> > ```
> > {: .bash}
> > 
> > ```
> > act wa
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


Note how the newly created directories and files are persistent, therefore can be re-accessed and re-used in future runs, even by containers instantiated from different images. All we have to do is to mount the filesystem image `my_overlay`.


### Run a Trinity genome assembly from inside the container

The directory we are in, `04_trinity`, contains sample inputs for a genome assembly, coming from one of the examples in the [Trinity Github repo]().

TBC


### Ephemeral writable containers

In some situations, you might need your container to be writable not to store persistent output files, but just to write temporary service files.  
*E.g.* this can happen with applications that want to write a dot-file in your home, such as a Python package, or containerised Jupyter notebooks that need to write runtime information under `/run`.  
In this context, a persistent overlay filesystem seems an overshooting. There alternative, simpler ways to set this up.

Singularity has a flag for rendering containers from SIF image files ephemerally writable. `--writable-tmpfs` will allocate a small amount of RAM for this purpose (configured by the sys admins, by default just a bunch of MB), e.g.:

```
$ singularity exec --writable-tmpfs library://ubuntu:18.04 touch ~/write-to-home
```
{: .bash}

Unless `$HOME` is bind mounted to the container (for security reasons it shouldn't), the newly created file will be gone after the container exits:

```
$ ls ~/write-to-home
```
{: .bash}

```
ls: /home/ubuntu/write-to-home: No such file or directory
```
{:. output}

There are situations where `--writable-tmpfs` is not usable, in particular if you are trying to write to a directory owned by *root*, such as `/run`.  
In this case, the solution is to create a host directory and bind mount it as the path you need to write into, *e.g.*:

```
$ mkdir ~/my_run
$ singularity exec -B ~/my_run:/run library://ubuntu:18.04 touch /run/running-file
```
{: .bash}

In this case, the file will also persist in the host directory after the container exits.
