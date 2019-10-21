---
title: "Sharing files with the host with Docker"
teaching: 15
exercises: 15
questions:
objectives:
- Learn how to mount host directories in a container
- Run a real-world bioinformatics application in a Docker container
keypoints:
- "Map host directories in the containers with the flag `-v <host dir>:<container dir>`"
- "Change working directory at container runtime with the flag `-w <working dir>`"
---

### Default directory access in Docker ###

Try and run the following to get to know what is the starting point in the Ubuntu container and what it contains. First: 

```
$ docker run ubuntu pwd
```
{: .bash}

```
/
```
{: .output}

Then:

```
docker run ubuntu ls -l
```
{: .bash}

```
total 64
drwxr-xr-x   2 root root 4096 Nov 12 20:56 bin
drwxr-xr-x   2 root root 4096 Apr 24  2018 boot
drwxr-xr-x   5 root root  340 Dec 19 08:01 dev
drwxr-xr-x   1 root root 4096 Dec 19 08:01 etc
drwxr-xr-x   2 root root 4096 Apr 24  2018 home
drwxr-xr-x   8 root root 4096 Nov 12 20:54 lib
drwxr-xr-x   2 root root 4096 Nov 12 20:55 lib64
drwxr-xr-x   2 root root 4096 Nov 12 20:54 media
drwxr-xr-x   2 root root 4096 Nov 12 20:54 mnt
drwxr-xr-x   2 root root 4096 Nov 12 20:54 opt
dr-xr-xr-x 125 root root    0 Dec 19 08:01 proc
drwx------   2 root root 4096 Nov 12 20:56 root
drwxr-xr-x   1 root root 4096 Nov 19 21:20 run
drwxr-xr-x   1 root root 4096 Nov 19 21:20 sbin
drwxr-xr-x   2 root root 4096 Nov 12 20:54 srv
dr-xr-xr-x  13 root root    0 Dec 14 13:27 sys
drwxrwxrwt   2 root root 4096 Nov 12 20:56 tmp
drwxr-xr-x   1 root root 4096 Nov 12 20:54 usr
drwxr-xr-x   1 root root 4096 Nov 12 20:56 var
```
{: .output}

You are in the root `/` directory of the container, and if you compare the listing of directories with what you get in the host (type `ls -l /` for this), you will notice the two are different; even directories with the same name will have in general different timestamps, suggesting they are in fact distinct directories.

This is the consequence of a default Docker behaviour: a container hasn't got any access to directories in the host filesystem (i.e. directories in the computer where you're running the container from)


### Accessing host directories ###

Docker has the ability to mount host directories into a container.  This allows you to add data to your container, as well as specify output directories you can use to store data after a container ends.  This is extremely useful as it's a bad idea to package up your containers with lots of data; it increases the size of the containers and makes them less portable (what if someone else wants to run the same container with different data?).

![Docker Volumes]({{ page.root }}/fig/docker-volume.png)

The docker daemon has a parameter called volume (`-v` or `--volume`), which we'll use to specify directories to be mounted.

The format is `-v /host/path:/container/path`.  Docker will create the directory inside the container if required.  Be aware the behaviour is different if you use absolute or relative paths, we use absolute paths here.

As an example, let us run the following:

```
$ docker run -v `pwd`:/data ubuntu ls -l /data
```
{: .bash}

Here we are using ``` `pwd` ``` as a shortcut for the current working directory. As a result of using the mapping option `-v`, the `ls` command run inside the container will display the content of the current directory in the host.

The `-v` flag maps host directories in the container, allowing to read/write within them. Let us use a container to create a file in a mapped directory:

```
$ docker run -v `pwd`:/data ubuntu touch /data/container1
```
{: .bash}

Now, let us look for that file in the host:

```
$ ls -l container1 
```
{: .bash}

```
-rw-r--r-- 1 root root 0 Dec 19 08:16 container1
```
{: .output}

The file created in the container is actually available from the host, as a consequence of volume mapping.


### Changing working directory at runtime

Docker has a flag to change working directory in the container, `-w` or `--workdir`. For instance let us use it to change dir to the mapped host directory, which allows us to adopt relative rather than absolute paths for files:

```
$ docker run -v `pwd`:/data -w /data ubuntu touch container2
$ ls -l container2
```
{: .bash}

```
-rw-r--r-- 1 root root 0 Dec 19 08:19 container2
```
{: .output}

This can be useful to make your workflow uniform, as different container providers may have different default working directories.


### More on volumes ###

Docker has several ways to mount data into containers. Here we've only partially covered the first one:

1. **bind mounts**: map a host directory inside the container. There are two possible syntaxes for this option, `-v` (or `--volume`) and `--mount`, the most significant difference being that `-v` is able to create the host directory to be mapped, if this doesn't exist, whereas `--mount` will throw an error. Docker is currently promoting `--mount` as the preferred syntax for mounting data.

2. **Docker volumes**: use storage spaces completely managed by Docker; they offer extra features compared to bind mounts.

3. **tmpfs mounts**: store data temporarily in the host memory.

[Manage data in Docker](https://docs.docker.com/storage/) contains detailed information on these options.


> ## Running BLAST from a container with Docker ##
> 
> We'll be running a BLAST (Basic Local Alignment Search Tool) example with a container from [BioContainers](https://biocontainers.pro).  BLAST is a tool bioinformaticians use to compare a sample genetic sequence to a database of known sequences; it's one of the most widely used bioinformatics tools.
> 
> To begin, try and pull the BLAST image `biocontainers/blast:v2.2.31_cv2` (this will take a little bit):
> 
> > ## Solution ##
> > 
> > ```
> > $ docker pull biocontainers/blast:v2.2.31_cv2
> > ```
> > {: .bash}
> > 
> > ```
> > v2.2.31_cv2: Pulling from biocontainers/blast
> > 22dc81ace0ea: Pull complete 
> > 1a8b3c87dba3: Pull complete 
> > [..]
> > Digest: sha256:238717ec69830ec62a19fc05c6f70183f218a13f7678864060f0157dc63dc54f
> > Status: Downloaded newer image for biocontainers/blast:v2.2.31_cv2
> > ```
> > {: .output}
> {: .solution}
> 
> Now, run a simple command to verify the image works, for instance `blastp -help`:
> 
> > ## Solution ##
> > 
> > ```
> > $ docker run biocontainers/blast:v2.2.31_cv2 blastp -help
> > ```
> > {: .bash}
> > 
> > ```
> > USAGE
> >   blastp [-h] [-help] [-import_search_strategy filename]
> > [..]
> >  -use_sw_tback
> >    Compute locally optimal Smith-Waterman alignments?
> > ```
> > {: .output}
> {: .solution}
> 
> The `03_blast` demo directory contains a human prion FASTA sequence, `P04156.fasta`, 
> and a gzipped reference database to blast against, `zebrafish.1.protein.faa.gz`. 
> Let us `cd` to it, and uncompress the database:
> 
> ```
> $ cd <top-level>/demos/03_blast
> $ gunzip zebrafish.1.protein.faa.gz
> ```
> {: .bash}
> 
> We need to prepare the zebrafish database with `makeblastdb` for the search, using the following command through a container:
> 
> ```
> $ makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> ```
> {: .bash}
> 
> To run it via Docker, you will need to mount the current directory to be able to read inputs and write outputs on the host.
> 
> Hint: the default directory in this BLAST image is `/data`, which can then be used as a convenient mount directory.
> 
> > ## Solution ##
> > 
> > ```
> > $ docker run -v `pwd`:/data biocontainers/blast:v2.2.31_cv2 makeblastdb -in zebrafish.1.protein.faa -dbtype prot
> > ```
> > {: .bash}
> {: .solution}
> 
> After the container has terminated, you should see several new files in the current directory.  We can now do the final alignment step using `blastp`, i.e.:
> 
> ```
> $ blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt
> ```
> {: .bash}
> 
> > ## Solution ##
> > 
> > ```
> > $ docker run -v `pwd`:/data biocontainers/blast:v2.2.31_cv2 blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt
> > ```
> > {: .bash}
> {: .solution}
> 
> The final results are stored in `results.txt`:
> 
> ```
> $ less results.txt
> ```
> {: .bash}
> 
> ```
>                                                                       Score     E
> Sequences producing significant alignments:                          (Bits)  Value
> 
>   XP_017207509.1 protein piccolo isoform X2 [Danio rerio]             43.9    2e-04
>   XP_017207511.1 mucin-16 isoform X4 [Danio rerio]                    43.9    2e-04
>   XP_021323434.1 protein piccolo isoform X5 [Danio rerio]             43.5    3e-04
>   XP_017207510.1 protein piccolo isoform X3 [Danio rerio]             43.5    3e-04
>   XP_021323433.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
>   XP_009291733.1 protein piccolo isoform X1 [Danio rerio]             43.5    3e-04
>   NP_001268391.1 chromodomain-helicase-DNA-binding protein 2 [Dan...  35.8    0.072
> [..]
> ```
> {: .output}
> 
> We can see that several proteins in the zebrafish genome match those in the human prion (interesting?).
{: .challenge}


### Useful container registries ###

There are a lot of applications (not just bioinformatics) already wrapped up in container images. 
Here's a small list of some of the web registries we use at Pawsey: 

* [Docker Hub](https://hub.docker.com)
* [Biocontainers](https://biocontainers.pro)
* [Quay](https://quay.io)^
* [Nvidia GPU Cloud (NGC)](https://ngc.nvidia.com)^

^The last two require you to create an account and login to access containers.


> ## Best practices ##
> 
> * Figuring out a standard way to consistently map host directories in container can help scripting and automation. For instance:
>     * ``` -v `pwd`:/data -w /data ``` can be useful when just working in the current directory
>     * ``` -v <DATA-DIRECTORY>:/data -w /data ``` can be useful if your workstation/cluster is organised with one directory called `<DATA-DIRECTORY>` that contains all sample data and reference data
> * Eventually, multiple volume mappings are allowed at the same time, for instance:
>         ```-v `pwd`:/data -v /reference-database:/ref -w /data ```
> * These syntaxes look ugly, but once learnt it can be reused with minimal variations
{: .callout}
