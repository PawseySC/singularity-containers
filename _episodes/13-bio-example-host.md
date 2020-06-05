---
title: "Share files with the host: BLAST, a bioinformatics demo"
teaching: 5
exercises: 15
questions:
objectives:
- Mount host directories in a container
- Pass specific variables to the container
- Run a real-world bioinformatics application in a container
keypoints:
- By default Singularity mounts the host current directory, and uses it as the container working directory
- Map additional host directories in the containers with the flag `-B`, or the variable SINGULARITY_BINDPATH
- Avoid mounting the `$HOME` directory, to better protect your sensitive data in the host
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


> ## And by the way, can we write inside a container?
> 
> Try and create a file called `example` in the container root directory.  (**Hint**: run `touch /example` inside the container).
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec docker://ubuntu:18.04 touch /example
> > ```
> > {: .bash}
> > 
> > ```
> > touch: cannot touch '/example': Read-only file system
> > ```
> > {: .output}
> > 
> > We have just learn something more on containers: by default, they are **read-only**.  How can we get a container to write files then?  Read on...
> {: .solution}
{: .challenge}


To summarise what we've learnt in the previous examples, we may say that a container ships an application and its dependencies by encapsulating them in an isolated, read-only filesystem.  In order for a container to access directories from the host filesystem (and write files), one needs to explicitly bind mount them.  The main exception here is the current work directory, which is bind mounted by default.


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
11-containers-intro.md  22-build-docker.md      33-gpu-gromacs.md       45-docker.md
12-singularity-intro.md 23-web-rstudio.md       41-workflow-engines.md  46-compose-web.md
13-bio-example-host.md  24-ml-python.md         42-x11-gnuplot.md       47-other-tools.md
14-build-intro.md       31-mpi-openfoam.md      43-wrappers.md
21-build-deffile.md     32-writable-trinity.md  44-setup-singularity.md
```
{: .output}

Also, we can write files in a host dir which has been bind mounted in the container:

```
$ singularity exec -B $TUTO docker://ubuntu:18.04 touch $TUTO/_episodes/example
$ singularity exec -B $TUTO docker://ubuntu:18.04 ls $TUTO/_episodes/example
```
{: .bash}

```
/home/ubuntu/singularity-containers/_episodes/example
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

We're going to use an image for the most recent BLAST version from the `quay.io` registry, *i.e.* `quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4`.  First, we'll pull the image.  This should take a few minutes (unless you had pulled the image in advance):

```
$ singularity pull docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4
```
{: .bash}


> ## Bonus: search for the BLAST image on an online registry
> 
> **If time allows**, you might want to give it a go with looking for the container image yourself.  
> Start with the assumption that most bioinformatics packages can be found within the *BioContainers* project (this is the repo/name you'll be looking for), and are hosted in both [Quay](https://quay.io) and [Biocontainers](https://biocontainers.pro).  
> These two registries contain the same images, they just offer a slightly different user interface.  At the time of writing, *Quay* has a cleaner and more readable interface compared to *BioContainers*; hopefully this will change in the future.  
> 
> > ## Solution: Quay
> > 
> > * Go to https://quay.io (NO registration required!);
> > * Locate the *Search* field on the top right of the page (you might need to widen the browser window), and type `blast`;
> > * We want an image from `biocontainers`, so look for `biocontainers/blast` and click on it;
> > * Click on the *Tags* icon on the left, and scroll the list of images to look for the highest Blast version (`2.9.0` at the time of writing); among the multiple tags for this version, identify the most recent one;
> > * At the time of writing, the resulting image will be `quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4`;
> > * You can click on the *Fetch* icon at the rightmost side of the record, select *Pull by Tag*, and then copy the full image name in your clipboard.
> {: .solution}
> 
> > ## Solution: BioContainers
> > 
> > * Go to https://biocontainers.pro;
> > * Click on the *Registry* button on the top of the page;
> > * In the new page, type `blast` in the search field;
> > * You will need to scroll a bit to find the proper `BLAST` entry; click on it;
> > * The list of images here is quite rich, with entries for *Docker*, *Singularity* and *Conda*; consider only the *Docker* entries, look for the highest Blast version (`2.9.0` at the time of writing) and, among the multiple tags for this version, identify the most recent one (*Hint*: sort by *Modified date*).  You might need to widen your window to read the full image names and tags, and on some smaller screens you won't be able to; **^Alternative download**
> > * At the time of writing, the resulting image will be `quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4`;
> > * You can click on the *Copy* icon just at the right of the image name field, to copy the full image name in your clipboard (you will need to get rid of *docker pull*).
> > 
> > **^Alternative download**: *Singularity* image
> > 
> > * Pick the highest version and latest tag from the list of *Singularity* entries;
> > * At the time of writing, the resulting image is again `quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4`;
> > * Click on the *Copy* icon just at the right of the image name field, to copy the full image name in your clipboard.  In this case, this is a *URL* to download the SIF image straight away: to achieve this, on your shell you will execute `wget <PASTE THE IMAGE NAME FROM CLIPBOARD>`.
> {: .solution}
{: .challenge}


> ## Run a test command
>
> Let us run a simple command using the image we just pulled, for instance `blastp -help`, to verify that it actually works.
>
> > ## Solution
> >
> > ```
> > $ singularity exec blast_2.9.0--pl526h3066fca_4.sif blastp -help
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
> > $ singularity exec ../blast/blast_2.9.0--pl526h3066fca_4.sif makeblastdb -in zebrafish.1.protein.faa -dbtype prot
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
> > $ singularity exec -B $TUTO/demos/blast_db blast_2.9.0--pl526h3066fca_4.sif blastp -query P04156.fasta -db $TUTO/demos/blast_db/zebrafish.1.protein.faa -out results.txt
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

When you're done, quit the view by hitting the `q` button.
