---
title: "Share files with the host: BLAST, a bioinformatics demo"
teaching: 5
exercises: 15
questions:
objectives:
- Learn how to mount host directories in a container
- Learn how to pass specific variables to the container
- Run a real-world bioinformatics application in a container
keypoints:
- By default Singularity mounts the host current directory, and uses it as the container working directory
- Map additional host directories in the containers with the flag `-B`, or the variable SINGULARITY_BINDPATH
- By default Singularity passes all host variables to the container
- Pass specific shell variables to containers by prefixing them with SINGULARITYENV_
---


### Access to directories in the host machine

Let's start and `cd` into the root demo directory:

```
$ cd $TUTO/demos
```
{: .bash}

What directories can we access from the container?

First, let us assess what the content of the root directory `/` looks like from outside *vs* inside the container, to highlight the fact that a container runs on his own filesystem:

```
$ ls /
```
{: .bash}

```
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  scratch  shared  srv  sys  tmp  usr  var
```
{: .output}


Now let's look at the root directory when we're in the container

```
$ singularity exec docker://ubuntu:18.04 ls /
```
{: .bash}

```
bin  boot  data  dev  environment  etc	home  lib  lib64  media  mnt  opt  proc  root  run  sbin  singularity  srv  sys  tmp  usr  var
```
{: .output}


> ## In which directory is the container running?
>
> For reference, let's check the host first:
>
> ```
> $ pwd
> ```
> {: .bash}
>
> ```
> /home/ubuntu/singularity-containers/demos
> ```
> {: .output}
>
> Now inspect the container.  (**Hint**: you need to run `pwd` in the container)
>
> > ## Solution
> >
> > ```
> > $ singularity exec docker://ubuntu:18.04 pwd
> > ```
> > {: .bash}
> >
> > ```
> > /home/ubuntu/singularity-containers/demos
> > ```
> > {: .output}
> >
> > Host and container working directories coincide!
> {: .solution}
{: .challenge}


> ## Can we see the content of the current directory inside the container?
>
> Hopefully yes ...
>
> > ## Solution
> >
> > ```
> > $ singularity exec docker://ubuntu:18.04 ls
> > ```
> > {: .bash}
> >
> > ```
> > blast  blast_db  gromacs  lolcow  lolcow_docker  lolcow_hpccm  nextflow  openfoam  pull_big_images.sh  python  rstudio	singularity  trinity
> > ```
> > {: .output}
> >
> > Indeed we can!
> {: .solution}
{: .challenge}


> ## How about other directories in the host?
>
> For instance, let us inspect `$TUTO/_episodes`.
>
> > ## Solution
> >
> > ```
> > $ singularity exec docker://ubuntu:18.04 ls $TUTO/_episodes
> > ```
> > {: .bash}
> >
> > ```
> > ls: cannot access '/home/ubuntu/singularity-containers/_episodes': No such file or directory
> > ```
> > {: .output}
> >
> > Host directories external to the current directory are not visible!  How can we fix this?  Read on...
> {: .solution}
{: .challenge}


> ## What happens on Pawsey HPC systems?
> 
> This last example won't work as expected on Zeus, Magnus and other Pawsey HPC machines.  
> This is due to site defaults that are meant to make users' life easier. In particular, `/group` and `/scratch` get bind mounted by default.  
> If you want to experience this example on Pawsey HPC, you should first `unset SINGULARITY_BINDPATH`.
{: .callout}


### Bind mounting host directories

Singularity has the runtime flag `--bind`, `-B` in short, to mount host directories.

There is a long syntax, which allows to map the host dir onto a container dir with a different name/path, `-B hostdir:containerdir`.  
There is also a short syntax, that just mounts the dir using the same name and path: `-B hostdir`.

Let's use the latter syntax to mount `$TUTO` into the container and re-run `ls`.

```
$ singularity exec -B $TUTO docker://ubuntu:18.04 ls $TUTO/_episodes
```
{: .bash}

```
11-containers-intro.md    14-build-intro.md         23-web-rstudio.md         32-writable-containers.md 44-docker.md
12-singularity-intro.md   21-build-deffile.md       24-ml-python.md           33-gpu-gromacs.md         45-other-tools.md
13-bio-example-host.md    22-build-docker.md        31-mpi-openfoam.md        41-workflow-engines.md
```
{: .output}

Now we are talking!

If you need to mount multiple directories, you can either repeat the `-B` flag multiple times, or use a comma-separated list of paths, *i.e.*

```
-B dir1,dir2,dir3
```
{: .bash}

Also, if you want to keep the runtime command compact, you can equivalently specify directories to be bind mounted using the environment variable `SINGULARITY_BINDPATH`:

```
$ export SINGULARITY_BINDPATH="dir1,dir2,dir3"
```
{: .bash}


> ## Mounting $HOME
>
> Depending on the site configuration of Singularity, user home directories might or might not be mounted into containers by default.  
> We do recommend that you **avoid mounting home** whenever possible, to avoid sharing potentially sensitive data, such as SSH keys, with the container, especially if exposing it to the public through a web service.
>
> If you need to share data inside the container home, you might just mount that specific file/directory, *e.g.*
>
> ```
> -B $HOME/.local
> ```
> {: .bash}
>
> Or, if you want a full fledged home, you might define an alternative host directory to act as your container home, as in
>
> ```
> -B /path/to/fake/home:$HOME
> ```
> {: .bash}
>
> Finally, you should also **avoid running a container from your host home**, otherwise this will be bind mounted as it is the current working directory.
{: .callout}


### How about sharing environment variables with the host?

By default, shell variables are inherited in the container from the host:

```
$ export HELLO=world
$ singularity exec docker://ubuntu:18.04 bash -c 'echo $HELLO'
```
{: .bash}

```
world
```
{: .output}

There might be situations where you want to isolate the shell environment of the container; to this end you can use the flag `-C`, or `--containall`:  
(Note that this will also isolate system directories such as `/tmp`, `/dev` and `/run`)

```
$ export HELLO=world
$ singularity exec -C docker://ubuntu:18.04 bash -c 'echo $HELLO'
```
{: .bash}

```

```
{: .output}

If you need to pass only specific variables to the container, that might or might not be defined in the host, you can define variables that start with `SINGULARITYENV_`; this prefix will be automatically trimmed in the container:

```
$ export SINGULARITYENV_CIAO=mondo
$ singularity exec -C docker://ubuntu:18.04 bash -c 'echo $CIAO'
```
{: .bash}

```
mondo
```
{: .output}


### Running BLAST from a container

We'll be running a BLAST (Basic Local Alignment Search Tool) example with a container from [BioContainers](https://biocontainers.pro).  BLAST is a tool bioinformaticians use to compare a sample genetic sequence to a database of known sequences; it's one of the most widely used bioinformatics packages.  
This example is adapted from the [BioContainers documentation](http://biocontainers-edu.biocontainers.pro/en/latest/running_example.html).

We're going to use the BLAST image `biocontainers/blast:v2.2.31_cv2`.  First, we'll pull the image.  This might take up to about 10 minutes (unless you had pulled the image in advance):

```
$ singularity pull docker://biocontainers/blast:v2.2.31_cv2
```
{: .bash}


> ## Run a test command
>
> Let us run a simple command using the image we just pulled, for instance `blastp -help`, to verify that it actually works.
>
> > ## Solution
> >
> > ```
> > $ singularity exec blast_v2.2.31_cv2.sif blastp -help
> > ```
> > {: .bash}
> >
> > ```
> > USAGE
> >   blastp [-h] [-help] [-import_search_strategy filename]
> >
> > [..]
> >
> >  -use_sw_tback
> >    Compute locally optimal Smith-Waterman alignments?
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


Now, the demo directory `demos/blast` contains a human prion FASTA sequence, `P04156.fasta`, whereas another directory, `demos/blast_db`, contains a gzipped reference database to blast against, `zebrafish.1.protein.faa.gz`.  Let us `cd` to the latter directory and uncompress the database:

```
$ cd $TUTO/demos/blast_db
$ gunzip zebrafish.1.protein.faa.gz
```
{: .bash}


> ## Prepare the database
>
> We then need to prepare the zebrafish database with `makeblastdb` for the search, using the following command through a container:
>
> ```
> $ makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> ```
> {: .bash}
>
> Try and run it via Singularity.
>
> > ## Solution
> >
> > ```
> > $ singularity exec blast_v2.2.31_cv2.sif makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> > ```
> > {: .bash}
> > ```
> > Building a new DB, current time: 11/16/2019 19:14:43
> > New DB name:   /home/ubuntu/singularity-containers/demos/blast_db/zebrafish.1.protein.faa
> > New DB title:  zebrafish.1.protein.faa
> > Sequence type: Protein
> > Keep Linkouts: T
> > Keep MBits: T
> > Maximum file size: 1000000000B
> > Adding sequences from FASTA; added 52951 sequences in 1.34541 seconds.
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


After the container has terminated, you should see several new files in the current directory (try `ls`).  
Now let's proceed to the final alignment step using `blastp`. We need to cd into `demos/blast`:

```
$ cd ../blast
```
{: .bash}


> ## Run the alignment
>
> Adapt the following command to run into the container:
>
> ```
> $ blastp -query P04156.fasta -db $TUTO/demos/blast_db/zebrafish.1.protein.faa -out results.txt
> ```
> {: .bash}
>
> Note how we put the database files in a separate directory on purpose, so that you will need to bind mount its path with Singularity.  Give it a go with building the syntax to run the `blastp` command.
>
> > ## Solution
> >
> > ```
> > $ singularity exec -B $TUTO/demos/blast_db blast_v2.2.31_cv2.sif blastp -query P04156.fasta -db $TUTO/demos/blast_db/zebrafish.1.protein.faa -out results.txt
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

The final results are stored in `results.txt`:

```
$ less results.txt
```
{: .bash}

```
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  XP_017207509.1 protein piccolo isoform X2 [Danio rerio]             43.9    2e-04
  XP_017207511.1 mucin-16 isoform X4 [Danio rerio]                    43.9    2e-04
  XP_021323434.1 protein piccolo isoform X5 [Danio rerio]             43.5    3e-04
  XP_017207510.1 protein piccolo isoform X3 [Danio rerio]             43.5    3e-04
  XP_021323433.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
  XP_009291733.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
  NP_001268391.1 chromodomain-helicase-DNA-binding protein 2 [Dan...  35.8    0.072
[..]
```
{: .output}
