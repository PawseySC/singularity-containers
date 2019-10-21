---
title: "Bioinformatics meets RStudio in containers"
teaching: 0
exercises: 20
questions:
objectives:
- Deploy a customised RStudio container for bioinformatics
keypoints:
- Containers are great way to manage R workflows.  You likely still want to have a local installation of R/Rstudio for some testing, but if you have set workflows, you can use containers to manage them.  You can also provide Rstudio servers for collaborators
- Also, docker-compose is a great way to manage complex Docker commands, as well as coordinating multiple containers
---

### RStudio example ###

R is a popular language in several domains of science, particularly because of its statistical packages.  It often requires installing a large number of dependencies, and installing these on an HPC system can be tedious.

Instead we can use an R container to simplify the process.


### Rocker ###

The group [Rocker](https://hub.docker.com/r/rocker) has published a large number of R images we can use, including an Rstudio image.  To begin, we'll pull a Tidyverse container image (contains R, RStudio, data science packages):

```
$ docker pull rocker/tidyverse:3.5
```
{: .bash}

We can now start this up:

```
$ docker run -d -p 80:8787 --name rstudio -v `pwd`/data:/home/rstudio/data -e PASSWORD=<Pick your password> rocker/tidyverse:3.5
```
{: .bash}

Here we're opening up the container port `8787` and mapping it to the host port `80` so we can access the Rtudio server remotely. Note you need to store a password in a variable; it will be required below for the web login.

You just need to open a web browser and point it to `localhost` if you are running Docker on your machine, or `<Your VM's IP Address>` if you are running on a cloud service.

You should see a prompt for credentials, with user defaulting to `rstudio`, and password..

Once you're done, stop the container with:

```
$ docker stop rstudio
```
{: .bash}


### Using RStudio images ###

The above example only provides a bare-bones RStudio image, but now we want to actually use some R packages.  The following example is based on a bioinformatics workshop at [OzSingleCell2018](https://github.com/IMB-Computational-Genomics-Lab/SingleCells2018Workshop).  We'll use their data for our Docker/Rstudio example.

To begin, let's clone the data (A trimmed down repo with their data has been created for this tutorial)

```
$ git clone https://github.com/skjerven/rstudio_ex.git
$ cd rstudio_ex
```
{: .bash}

For this example, we'll use an RStudio image thas has already been built.  R images can take a while to build sometimes, depending on the number of packages and dependencies you're installing.  The Dockerfile used here is included, and we'll go through it to explain how Docker builds images.
 
```
FROM rocker/tidyverse:3.5

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
      autoconf \
      automake \
      g++ \
      gcc \
      gfortran \
      make \
      && apt-get clean all \
      && rm -rf /var/lib/apt/lists/*

RUN mkdir -p $HOME/.R
COPY Makevars /root/.R/Makevars

RUN Rscript -e "library('devtools')" \
      -e "install_github('Rdatatable/data.table', build_vignettes=FALSE)" \
      -e "install.packages('reshape2')" \
      -e "install.packages('fields')" \
      -e "install.packages('ggbeeswarm')" \
      -e "install.packages('gridExtra')" \
      -e "install.packages('dynamicTreeCut')" \
      -e "install.packages('DEoptimR')" \
      -e "install.packages('http://cran.r-project.org/src/contrib/Archive/robustbase/robustbase_0.90-2.tar.gz', repos=NULL, type='source')" \
      -e "install.packages('dendextend')" \
      -e "install.packages('RColorBrewer')" \
      -e "install.packages('locfit')" \
      -e "install.packages('KernSmooth')" \
      -e "install.packages('BiocManager')" \
      -e "source('http://bioconductor.org/biocLite.R')" \
      -e "biocLite('Biobase')" \
      -e "biocLite('BioGenerics')" \
      -e "biocLite('BiocParallel')" \
      -e "biocLite('SingleCellExperiment')" \
      -e "biocLite('GenomeInfoDb')" \
      -e "biocLite('GenomeInfgoDbData')" \
      -e "biocLite('DESeq')" \
      -e "biocLite('DESeq2')" \
      -e "BiocManager::install(c('scater', 'scran'))" \
      -e "library('devtools')" \
      -e "install_github('IMB-Computational-Genomics-Lab/ascend', ref = 'devel')" \
      && rm -rf /tmp/downloaded_packages
```
{: .source}

The first line, `FROM`, specifies a base image to use.  We could build up a full R image from scratch, but why waste the time.  We can use Rocker's pre-built image to start with and simplify our lives.

`RUN apt-get update` is installing some packages we'll need via Ubuntu's package manager.  Really all we're installing here are compilers.

The next section adds some flags and options we want to use when building R packages, by copying a file from the build context, `Makvevars`.

The last section is the main R package installation section.  Here we run several different installation methods:

* `install.packages()` is R's standard method for installing packages from CRAN.  We also install the `robustbase` package from source here.
* `biocLite()` is [Bioconductor's](https://bioconductor.org) method for installing packages
* `BiocManager` is a [CRAN package](https://cran.r-project.org/package=BiocManager) for installing bioinformatics software
* `install_github()` is method for installing R packages from GitHub.

We'll skip building this image for now, and just pull and use a prebuilt image.  We're also going to use `docker-compose` to help with setting up our container (see previous episode on long running servies). Here we'll use it for managing several options we want to use for our Rstudio image.

```
version: "2"

services:
  rstudio:
    restart: always
    image: bskjerven/oz_sc:latest
    container_name: rstudio
    volumes:
      - "$HOME/rstudio_ex/data:/home/rstudio/data"
    ports:
      - 80:8787
    environment:
      - USER=rstudio
      - PASSWORD=rstudiopwd
```
{: .source}

This yaml file simple tells Docker what image we want to run along with some options (like which volumes to mount, username/password, and what network ports to use).

To begin, make sure you're in the `rstudio_ex` directory in your home (where we cloned the repo).  Simply type:

```
$ docker-compose up
```
{: .bash}

Docker will pull the `oz_sc:latest` image first (if it's not present on your system yet); once that's complete you'll see output from the RStudio server:

```
[..]
Recreating rstudio ... done
Attaching to rstudio
rstudio    | [fix-attrs.d] applying owners & permissions fixes...
rstudio    | [fix-attrs.d] 00-runscripts: applying...
rstudio    | [fix-attrs.d] 00-runscripts: exited 0.
rstudio    | [fix-attrs.d] done.
rstudio    | [cont-init.d] executing container initialization scripts...
rstudio    | [cont-init.d] add: executing...
rstudio    | Nothing additional to add
rstudio    | [cont-init.d] add: exited 0.
rstudio    | [cont-init.d] userconf: executing...
rstudio    | [cont-init.d] userconf: exited 0.
rstudio    | [cont-init.d] done.
rstudio    | [services.d] starting services
rstudio    | [services.d] done.
```
{: .output}

This is annoying, though...we need our terminal back.  Luckily, Docker lets you run processes in the background.  Kill the RStudio process with `CTRL-C`, and the rerun `docker-compose` with the `-d` flag:

```
$ docker-compose up -d
```
{: .bash}

Shortly after that starts, open a web browser and go to `localhost` if you are running Docker on your machine, or `<Your VM's IP Address>` if you are running on a cloud service. You should see an Rstudio login, and we've set the username to `rstudio` and password to `rstudiopwd`.

Once logged in, you type (note this is the R shell):

```
> source('data/SC_script.r')
```
{: .r}

to run the tutorial (it may take a few minutes).  We can refer to the [OzSingleCell2018](https://github.com/IMB-Computational-Genomics-Lab/SingleCells2018Workshop) repo for details on each step.

To stop your Rstudio image, simply type from the `rstudio_ex` directory:

```
$ docker-compose down
```
{: .bash}


### Running a scripted R workflow on HPC with Shifter ###

We can run the same analysis on HPC through command line using Shifter. We can use the same container image, but rather than an RStudio GUI we'll use the `Rscript` command to execute the script.

To get started let's pull the required R container image:

```
$ module load shifter
$ sg $PAWSEY_PROJECT -c 'shifter pull bskjerven/oz_sc:latest'
```
{: .bash}

Now let's change directory to either `$MYSCRATCH` or `$MYGROUP`, e.g.

```
$ cd $MYSCRATCH
```
{: .bash}

With your favourite text editor, create a SLURM script, we'll call it `rscript-bio.sh` (remember to specify your Pawsey project ID in the script!):

```
#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --job-name=rstudio-bio

module load shifter

# clone Git repo with sample data and script
git clone https://github.com/skjerven/rstudio_ex.git
cd rstudio_ex

# run R script
srun --export=all shifter run bskjerven/oz_sc:latest Rscript data/SC_Rscript.r
```
{: .bash}

Let's submit the script via SLURM:

```
$ sbatch --reservation <your-pawsey-reservation> rscript-bio.sh
```
{: .bash}
