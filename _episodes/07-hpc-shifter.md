---
title: "Run containers on HPC with Shifter (and Singularity)"
teaching: 15
exercises: 5
questions:
objectives:
- Learn how to manage and run containers on a HPC cluster with Shifter
keypoints:
- "Shifter has a quite simple syntax that allows to pull, manage and run containers on HPC systems"
- "`shifter pull` and `shifter run` are the key commands"
---

### Why not Docker on HPC? ###

There are a few issues preventing Docker from being used as a container engine on HPC systems:

* Security: Docker requires root privileges
* Batch systems: doesn't integrate well with schedulers
* Underlying kernel: usually requires an up-to-date kernel

Fortunately, a number of alternatives are available to run containers at HPC facilities, including:

* [Shifter](https://docs.nersc.gov/programming/shifter/overview/): developed by NERSC and Cray, Docker-like interface, MPI support 
* [CSCS Shifter](https://user.cscs.ch/tools/containers/): forked by CSCS, adds features including GPU support 
* [Singularity](https://www.sylabs.io/singularity/): originally developed by LBL, has its own image format and can run Docker containers as well

At the moment, Pawsey is using CSCS Shifter on its HPC systems, and therefore this will be the tool of choice in this tutorial.

NCI (the National Computational Infrastructure in Canberra) is using Singularity on its HPC systems, for which examples will be provided as well.


> ## How to login on Pawsey HPC systems? ##
> 
> Pawsey currently has two sytems, **Magnus** and **Zeus**. We're using Zeus for this tutorial. You can login using the `ssh` command and your Pawsey access credentials (they will be provided for live workshops):
> 
> ```
> $ ssh <your-pawsey-account-name>@zeus.pawsey.org.au
> ```
> {: .bash}
> 
> After this, you'll be prompted to enter your Pawsey account password.
{: .callout}


### Pulling and managing images with Shifter ###

To use Shifter on Pawsey HPC systems, we need first to load the corresponding module:

```
$ module load shifter
```
{: .bash}

In principle, the command to pull container images is very similar to Docker, `shifter pull`.  
However, to avoid disk quota issues on Pawsey HPC systems, the following syntax is recommended, that makes use of the `sg` linux command, for instance:

```
$ sg $PAWSEY_PROJECT -c 'shifter pull ubuntu'
```
{: .bash}

```
# image     : index.docker.io/library/ubuntu/latest
# cacheDir  : /group/shifterrepos/mdelapierre/.shifter/cache
# tmpDir    : /tmp
# imageDir  : /group/shifterrepos/mdelapierre/.shifter/images
> save image layers ...
> pulling        : sha256:f85999a86bef2603a9e9a4fa488a7c1f82e471cbb76c3b5068e54e1a9320964a
> pulling        : sha256:da1315cffa03c17988ae5c66f56d5f50517652a622afc1611a8bdd6c00b1fde3
[..]
> extracting     : /group/shifterrepos/mdelapierre/.shifter/cache/sha256:f85999a86bef2603a9e9a4fa488a7c1f82e471cbb76c3b5068e54e1a9320964a.tar
> make squashfs ...
> create metadata ...
# created: /group/shifterrepos/mdelapierre/.shifter/images/index.docker.io/library/ubuntu/latest.squashfs
# created: /group/shifterrepos/mdelapierre/.shifter/images/index.docker.io/library/ubuntu/latest.meta
```
{: .output}

```
$ sg $PAWSEY_PROJECT -c 'shifter pull centos'
$ sg $PAWSEY_PROJECT -c 'shifter pull busybox'
```
{: .bash}

Similar again to Docker, we can list locally pulled images with `shifter images`:

```
$ shifter images
```
{: .bash}

```
library/centos                   latest                       ea4b646d9000   2018-11-27T07:05:23   69.62MB      index.docker.io
library/busybox                  latest                       7dc9d60af829   2018-12-19T22:31:48   704.00KB     index.docker.io
library/ubuntu                   latest                       d71fc6939e16   2018-12-19T22:30:41   29.94MB      index.docker.io
```
{: .output}

and remove undesired images with `shifter rmi`:

```
$ shifter rmi busybox
```
{: .bash}

```
removed image index.docker.io/library/busybox/latest
```
{: .output}


### Running images with Shifter ###

Let us change directory to our group directory with:

```
$ cd $MYGROUP
```
{: .bash}

and then run `ls` using the Ubuntu image we just pulled, via `shifter run`:

```
$ shifter run ubuntu ls
```
{: .bash}

The output will display the content of the current host directory!

A few differences in behaviour can be noticed compared to Docker, such that using Shifter typically requires to specify less options and flags:

* by default, some relevant directories in the Pawsey HPC filesystems are mounted in the containers; these include `/group`, `/scratch`, `/pawsey` and `/tmp`.  
  **Note**: `/home` is NOT mounted instead;
* if running from a mapped host directory, this becomes the working directory at container runtime;
* standard input is always open, allowing redirection;
* the host user is automatically set for the container;
* Shifter automatically removes containers after execution is terminated.

As additional examples, you might want to run:

```
$ shifter run ubuntu ls /
$ shifter run ubuntu whoami
```
{: .bash}

Note how no flag is required to run a container interactively. To launch an interactive shell from within the container, just run it without any commands, for instance:

```
$ shifter run ubuntu
```
{: .bash}

```
mdelapierre@zeus-1:/group/pawsey0001/mdelapierre$ exit   # or hit CTRL-D
```
{: .bash}

Shifter has support to run containers exploiting MPI parallelism and GPU acceleration (the latter only through CSCS Shifter).

* `shifter run --mpi` allows containers to take advantage on inter-node communication on the host fabric. The container image needs to have been built with MPI libraries thatare ABI compatible with the host MPI libraries;

* to run GPU enabled containers, no extra flags are required.


### Using Shifter with a job scheduler ###

Shifter is compatible with **SLURM**, the job scheduler installed on Pawsey HPC systems. In particular, SLURM job executor `srun` is compatible with `shifter run`, and the two syntaxes can be combined together.

As an example, the following script uses a Ubuntu container to output the machine hostname:

```
#!/bin/bash -l
  
#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --export=NONE
#SBATCH --job-name=blast

module load shifter

srun --export=all shifter run ubuntu hostname
```
{: .bash}

Now use your favourite text editor to copy paste the script above in a file called `hostname.sh` somewhere under `$MYSCRATCH` or `$MYGROUP` (remember to specify your Pawsey Project ID in the script!),

and then submit this script using SLURM. If you are running this during a live workshop, use the flag `--reservation <your-pawsey-reservation>` to use the compute nodes that have been reserved for the event:

```
$ sbatch --reservation <your-pawsey-reservation> hostname.sh
```
{: .bash}


### Can Shifter build container images? ###

Shifter does not allow to build container images. The best way to create an image to be pulled and run through it is to use Docker on a distinct machine (see previous episode).


> ## Run a Python app in a container on HPC ##
> 
> First, pull the container `continuumio/miniconda3:4.5.12`.
> 
> Then, with your favourite text editor create a file called `app.py` with the following content:
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
> Now, run `python app.py` using the the container image you have just pulled. For instance, give the input filename as an argument to the app.
> 
> Finally, re-run it by means of a SLURM script called `python_slurm.sh`.
> 
> > ## Solution ##
> > 
> > Pull the container image:
> > 
> > ```
> > $ sg $PAWSEY_PROJECT -c 'shifter pull continuumio/miniconda3:4.5.12'
> > ```
> > {: .bash}
> > 
> > Run the app:
> > 
> > ```
> > $ shifter run continuumio/miniconda3:4.5.12 python app.py input
> > ```
> > {: .bash}
> > 
> > SLURM script for scheduler submission, `python_slurm.sh` (insert Pawsey Project ID!):
> > 
> > ```
> > #!/bin/bash -l
> > 
> > #SBATCH --account=<your-pawsey-project>
> > #SBATCH --partition=workq
> > #SBATCH --ntasks=1
> > #SBATCH --time=00:05:00
> > #SBATCH --export=NONE
> > #SBATCH --job-name=python
> > 
> > module load shifter
> > 
> > srun --export=all shifter run continuumio/miniconda3:4.5.12 python app.py input
> > ```
> > {: .bash}
> > 
> > SLURM submission:
> > 
> > ```
> > $ sbatch --reservation <your-pawsey-reservation> python_slurm.sh
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


### HPC containers with Singularity at NCI ###

Raijin, the NCI HPC system, uses Singularity as the container engine, instead of Shifter. However, much of the user-facing interface is similar, or even the same. The biggest difference is that on Raijin you cannot pull and use your own images; instead, you’ll need to contact the NCI Help Desk at <mailto:help@nci.org.au> and ask for your image to be added to the library.
 
This allows the NCI staff to inspect and sanitise the containers before use. For example, to ensure that they allow the use of system MPI libraries, or at least contain a compatible version. The images must be able to be built using a build script, rather than being distributed as just an opaque filesystem image. Mount points for the systems Lustre filesystems will also be included on build so that all of your user data is available at the same location as in the native image (e.g. `/home`, `/short`, and `/g/data`).
 
First of all, let us load the Singularity module:

```
module load singularity
```
{: .bash}

You can run an interactive shell inside the container using the `singularity shell` command. Here, we are using a CentOS image:
 
```
$ singularity shell /apps/singularity/images/centos7/centos7-latest.simg 
```
{: .bash}

```
Singularity: Invoking an interactive shell within container...
 
Singularity centos7-2018051701.simg:~> cat /etc/centos-release
CentOS Linux release 7.5.1804 (Core) 
Singularity centos7-2018051701.simg:~> 
```
{: .output}

```
Singularity centos7-2018051701.simg:~> whoami
```
{: .bash}

```
bjm900
```
{: .output}

```
Singularity centos7-2018051701.simg:~> exit
```
{: .bash}

You can run a specific command within the container using the `singularity exec` command:

```
$ singularity exec /apps/singularity/images/centos7/centos7-latest.simg cat /etc/centos-release
```
{: .bash}

```
CentOS Linux release 7.5.1804 (Core) 
```
{: .output}

Note that the CentOS version in the container image is 7.5, whereas Raijin’s native OS is (currently) CentOS 6.10:

```
$ cat /etc/centos-release
```
{: .bash}

```
CentOS release 6.10 (Final)
```
{: .output}

On Raijin, Singularity is also integrated with the **PBS** batch scheduling system. This allows you to specify the image to use via the directive `#PBS -l image` in your job script (or similarly on the `qsub` command line):

```
#!/bin/bash
#PBS -P <your-nci-project>
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=00:05:00
#PBS -l image=centos7

module load singularity

cat /etc/centos-release
```
{: .bash}

Of course, you can also just use `singularity exec` within your job script as in the example above:

```
#!/bin/bash
#PBS -P <your-nci-project>
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=00:05:00

module load singularity

singularity exec /apps/singularity/images/centos7/centos7-latest.simg cat /etc/centos-release
```
{: .bash}

