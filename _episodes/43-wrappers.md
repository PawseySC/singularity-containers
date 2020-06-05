---
title: "Streamline the user experience: bash wrappers and modules"
teaching: 10
exercises: 10
questions:
objectives:
- Simplify containers usage by means of bash wrappers
- Discuss how to deploy containers and their wrappers using modules
keypoints:
- It is possible to devise a quite general wrapper template for containerised application
- The key information to setup the wrappers is the container image, and the commands one needs to run from that image
- It is possible to write a minimal modulefile, that allows to setup the shell environment to use containerised applications through wrappers
---


### Can we standardise the use of containers, to simplify the required syntax?

To answer this question, let's grab the BLAST container we used in the demo on sharing files with the host (the image will be quickly pulled from cache if you ran that demo):

```
$ cd $TUTO/demos/wrap_blast
$ singularity pull docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4
```
{: .bash}

Now, let's think about the typical usage of a containerised application.  Once the container image is available in the local disk, in the vast majority of cases you'll use it to execute some command in this way:

```
singularity exec ./blast_2.9.0--pl526h3066fca_4.sif <CMD> <ARGS>
```
{: .source}

As a plain, useful example, let's suppose we want to get the help output from the `blastp` command:

```
$ singularity exec ./blast_2.9.0--pl526h3066fca_4.sif blastp -help
```
{: .bash}

We can break this into logical parts; let's write a script called `blastp.1` for convenience:

```
#!/bin/bash

image_dir="."
image_name="blast_2.9.0--pl526h3066fca_4.sif"

cmd="blastp"

args="$@"

singularity exec $image_dir/$image_name $cmd $args
```
{: .source}

Look at how general the expression in the last line of this script is!  
We're also using shell variables to express tool- and command- specific information.  Of these, the image location `image_dir` and name `image_name` are set at the time we pull the image.  The command name, `cmd`, might change from command to command.  So, for instance, we might write a script for the command `makeblastdb` by only changing that line:

```
#!/bin/bash

image_dir="."
image_name="blast_2.9.0--pl526h3066fca_4.sif"

cmd="makeblastdb"

args="$@"

singularity exec $image_dir/$image_name $cmd $args
```
{: .source}

How about the value we assigned to the command arguments variable, `args`?  Well, that's bash syntax.  If you execute this script, bash will assign to `$@` the full list of arguments that you append to the script in the command line.  
To see a practical example, let's make the `blastp.1` script executable (using `chmod`) and run it with the `-help` argument:

```
$ chmod +x blastp.1
```
{: .bash}

```
$ ./blastp.1 -help
```
{: .bash}

```
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
[..]
 -use_sw_tback
   Compute locally optimal Smith-Waterman alignments?
```
{: .output}

From the output, you can see that the `blastp` command actually got the `-help` flag right, and this was thanks to the usage of `$@` in the script.

So to summarise this section, we've written a simple bash script that wraps around the Singularity `exec` approach, so that to run `blastp` from a container you simply type:

```
$ ./blastp.1 <ARGS>
```
{: .bash}

Why the `.1` extension?  Well, this is just because the story is not over...


### A (quite) general bash wrapper for containerised applications

In the first iteration of a bash wrapper for containerised commands, we need to provide 3 pieces of information in the script: image location, image name and command name.  Can we further simplify and generalise this?

**Yes**.  With a couple of extra bash commands and assumptions, we can make it so that the only required information will be the **container image name**.

First, let's get rid of the command name.  
Let's assume that we're calling the wrapper with the same name of the command we want it to execute.  Then, we're going to use the bash variable `$0`; used inside a script, it contains the full path of the script itself; we're also using the bash command `basename`, that extract a file or directory name out of its full path.  The `cmd` variable becomes:

```
cmd="$(basename $0)"
```
{: .source}

And now, let's generalise the image location.  
Let's assume that we're storing the wrappers in the same directory where the image is located.  Then, we can use the bash command `dirname` to extract the location of a file or directory out of its full path.  The `image_dir` variable becomes:

```
image_dir="$(dirname $0)"
```
{: .source}

So we can now have a general bash wrapper for BLAST commands from the container image `blast_2.9.0--pl526h3066fca_4.sif`:

```
#!/bin/bash

image_dir="$(dirname $0)"
image_name="blast_2.9.0--pl526h3066fca_4.sif"

cmd="$(basename $0)"

args="$@"

singularity exec $image_dir/$image_name $cmd $args
```
{: .source}

To create a wrapper for `blastp`, all we have to do is to create a script named `blastp` with that content.  Then, we can do the same for `makeblastdb`, `blastn`, `blastx` and so on.  
To limit the number of files, we might even just have a single copy of this script, *e.g.* named `blastp`, and then create symbolic links for the other commands, for instance:

```
$ ln -s blastp makeblastdb
```
{: .bash}

What if we need bash wrappers for the Trinity assembler from the pulled image ?  
Well, just make a new script with a different `image_name`, named according to the required command:

```
image_name="trinityrnaseq_2.8.6.sif"
```
{: .source}


### How general is this approach?

Well, quite general probably.  It can be used every time you would use containers with this Singularity syntax:

```
singularity exec <IMAGE> <CMD> <ARGS>
```
{: .source}

This will also work with MPI containers and Slurm, as the corresponding syntax does not impact such form:

```
mpirun -n <NNODES> singularity exec <IMAGE> <CMD> <ARGS>
srun -n <NNODES> singularity exec <IMAGE> <CMD> <ARGS>
```
{: .source}

Of course there are some corner cases.  
For instance, for GPU enabled containers, after `exec` in the wrapper you will need to add `--nv`.  
Using overlays requires adding `--overlay <OVERLAY FILEPATH>`, with the file path possibly specified using a shell variable that you can define prior to executing the wrapper.  
Wrappers to launch GUI sessions will also require some tweaking.  


### What if we need to bind mount some host directories?

This is a case worth commenting in this context.  

Specifying the paths to be bind mounted as additional flags in the wrappers is not really general nor portable.  

So what you want to do here is to use `$SINGULARITY_BINDPATH`, defining the required paths prior to execution of the application.  
If you have a standard setup on your system, where all the data go under the same parent directory (*e.g.* `/data`), you might even want to define the variable in the startup scripts (`~/.bashrc`,...).  This can be quite a good practice in simplifying your production environment, and making it more robust.  
In this respect, in Pawsey HPC systems the singularity module adds `/group` and `/scratch` to the the bind path, so you don't have to worry about bind mounting data directories at all.


### Bonus: use modules to handle bash wrappers

So far in this episode, we've devised a scenario to deploy a containerised application in a streamlined way:

1. define the container image you need;
2. pull it in a directory;
3. in that same directory, create bash wrappers for the commands you need to execute from that container.

If you're in a system with lots of other applications, you might want to tidy up the environment by using modules.  Here, we're using the [Environment Modules](http://modules.sourceforge.net) implementation.  This tutorial provides a template installation [script]({{ page.root }}/files/install-modules.sh) for a Linux box.  
**Note that** discussing modules in details is off scope here, we're just using them to show how to organise containerised applications.

We have just said that all relevant files for our containerised application, *e.g.* BLAST, are in a single location.  
To run this example, there's already a directory made ready in your current work directory, `$TUTO/demos/wrap_blast`, namely `apps/blast/2.9.0/bin`.  It contains four bash wrappers:

```
$ ls apps/blast/2.9.0/bin
```
{: .bash}

```
blastn      blastp      blastx      makeblastdb
```
{: .output}

To get ready for this example, let us also pull the BLAST image there:

```
$ singularity pull --dir apps/blast/2.9.0/bin docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4
```
{: .bash}

Now, we can think of a minimal modulefile to setup BLAST in our environment:

```
#%Module1.0######################################################################
##
## blast modulefile
##
proc ModulesHelp { } {
    puts stderr "\tModule for blast version 2.9.0\n"
    puts stderr "\tThis module uses the container image docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4"
}

module-whatis   "edits the PATH to use the tool blast version 2.9.0"

prepend-path     PATH            $env(TUTO)/demos/wrap_blast/apps/blast/2.9.0/bin
```
{: .source}

In general, the string associated to `PATH` will need to be customised case by case, same as the `help` and `whatis` strings.  
A copy of this modulefile is under `modulefiles/` in the current path.  
Shall we try it?  Let's go!

Let's tell modules to look for modules in this directory, and then test it can find the blast module:

```
$ module use $(pwd)/modulefiles
$ module avail
```
{: .bash}

```
----------------------------------------------- /data/work/gitrepos/Trainings/singularity-containers/demos/wrap_blast/modulefiles ------------------------------------------------
blast/2.9.0  

------------------------------------------------------------------------- /usr/share/modules/modulefiles -------------------------------------------------------------------------
dot  module-git  module-info  modules  null  use.own  

```
{: .output}

It's there!  Let's `load` it then:

```
$ module load blast/2.9.0
```
{: .bash}

Can we now see the wrappers in there?

```
$ which blastp
```
{: .bash}

```
/home/ubuntu/singularity-containers/demos/wrap_blast/apps/blast/2.9.0/bin/blastp
```
{: .output}

Sure!  Let's test it with the usual `-help` flag:

```
$ blastp -help
```
{: .bash}

```
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
[..]
 -use_sw_tback
   Compute locally optimal Smith-Waterman alignments?
```
{: .output}

Containerised BLAST with wrappers and modules: the experience looks like a traditional installation!


### Final thoughts

So, we've shown you how to effectively hide containers under the hood to provide a simplified user experience, while gaining in reproducibility, portability, productivity and more.

Why bothering with learning the longer story of the Singularity syntax then?  Well, containers are a powerful technology, but also a complex one.  
Even if you're going to use them through a friendlier interface, it's still crucial to know how thing work underneath, to be aware of the corresponding limitations, and possibly also to be able to fix the setup when things go wrong.
