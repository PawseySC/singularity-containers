---
title: "Sharing files with the host with Docker"
teaching: 10
exercises: 5
questions:
objectives:
- Learn how to mount host directories in a container
keypoints:
- "Map host directories in the containers with the flag `-v <host dir>:<container dir>`"
---

### Directory and file defaults in Docker ###

Try and run the following to get to know what is the starting point in the Ubuntu container and what it contains:

```
$ docker run ubuntu pwd
```
{: .bash}

```
/
```
{: .output}

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

Now try and create an empty file and then see who is the owner (we're feeding two commands at once to the container by separating them with a semi-colon, and running through `bash -c`):

```
$ docker run ubuntu bash -c 'touch empty-file ; ls -l empty-file'
```
{: .bash}

```
-rw-r--r-- 1 root root 0 Dec 19 08:06 empty-file
```
{: .output}

The file is owned by the root user!

What we have just seen is a consequence of some Docker defaults:

* a container hasn't got any access to directories in the host filesystem (i.e. directories in the computer where you're running the container from)
* as by default a container is run as root, any created file is owned by the group user.


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

Finally, Docker has a flag to change working directory in the container, to avoid using full paths, `-w` or `--workdir`; for instance let us use it to change dir to the mapped host directory:

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


> ## Run a Python app in a container with I/O ##
> 
> With your favourite text editor create a file called `app.py` with the following content:
> 
> ```
> import sys
>   
> def print_sums(data):
>     with open("row_sums",'w') as output:
>         for line in data:
>             row = 0
>             for word in line.strip().split():
>                 row += int(word)
>             output.write(str(row)+"\n")
>             print("Sum of the row is ",row)
> 
> if len(sys.argv) > 1 and sys.argv[1] != "-":
>     with open(sys.argv[1], 'r') as infile:
>         print_sums(infile)
> else:
>     print_sums(sys.stdin)
> ```
> {: .python}
> 
> and an input file `input` containing:
> 
> ```
> 1 2 3
> 4 5 6
> 7 8 9
> ```
> {: .source}
> 
> The app reads rows containing integers and outputs their sums line by line. Input can be given through file or via standard input. The output is produced both in formatted form through standard output and in raw form written to a file named `row_sums`.
> 
> Now, run `python app.py` using the the container image `continuumio/miniconda3:4.5.12` you previously pulled. Give the input filename as an argument to the app.
> 
> > ## Solution ##
> > 
> > Run with input file as argument:
> > 
> > ```
> > $ docker run -v `pwd`:/data -w /data continuumio/miniconda3:4.5.12 python app.py input
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


> ## Best practices ##
> 
> * Figuring out a standard way to consistently map host directories in container can help scripting and automation. For instance:
>     * ``` -v `pwd`:/data -w /data ``` can be useful when just working in the current directory
>     * ``` -v /<DATA-DIRECTORY>:/data -w /data ``` can be useful if your workstation/cluster is organised with one directory called `<DATA-DIRECTORY>` that contains all sample data and reference data
> * Eventually, multiple volume mappings are allowed at the same time, for instance:
>         ```-v `pwd`:/data -v /reference-database:/ref -w /data ```
> * These syntaxes look ugly, but once learnt it can be reused with minimal variations
{: .callout}
