---
title: "RStudio deployment for fun and profit"
teaching: 0
exercises: 20
questions:
objectives:
- Run an R workflow both through RStudio and the terminal using containers
keypoints:
- Containers are great way to manage R workflows.  You likely still want to have a local installation of R/Rstudio for some testing, but if you have set workflows, you can use containers to manage them.  You can also provide Rstudio servers for collaborators
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


### Running a scripted R workflow on the shell ###

Let us create a dedicated directory for this example:

```
$ mkdir r_example
$ cd r_example
```
{: .bash}

We are going to use a minimalistic example taken from the workshop [Programming with R](http://swcarpentry.github.io/r-novice-inflammation/) by the Software Carpentry.
The script `readings-06.R` from their [Episode 5](http://swcarpentry.github.io/r-novice-inflammation/05-cmdline/index.html) is made available here for convenience, you can copy-paste the content in a file using your favourite text editor:

```
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- args[1]
  filenames <- args[-1]
  stopifnot(action %in% c("--min", "--mean", "--max"))

  if (length(filenames) == 0) {
    process(file("stdin"), action)
  } else {
    for (f in filenames) {
      process(f, action)
    }
  }
}

process <- function(filename, action) {
  dat <- read.csv(file = filename, header = FALSE)

  if (action == "--min") {
    values <- apply(dat, 1, min)
  } else if (action == "--mean") {
    values <- apply(dat, 1, mean)
  } else if (action == "--max") {
    values <- apply(dat, 1, max)
  }
  cat(values, sep = "\n")
}

main()
```
{: .r}

Let us download and unzip the required sample dataset:

```
$ wget http://swcarpentry.github.io/r-novice-inflammation/data/r-novice-inflammation-data.zip
$ unzip -q r-novice-inflammation-data.zip
```
{: .bash}

Now, we can run the R script using the R container we pulled; we're going to compute average values in this example:

```
$ docker run -v `pwd`:/data -w /data rocker/tidyverse:3.5 Rscript readings-06.R --mean data/inflammation-*.csv
```
{: .bash}


### Using an RStudio web server to run the analysis ###

Let us start up the web server using the following `docker` command:

```
$ docker run -d -p 80:8787 --name rstudio -v `pwd`/data:/home/rstudio/data -e PASSWORD=<Pick your password> rocker/tidyverse:3.5
```
{: .bash}

Here we're opening up the container port `8787` and mapping it to the host port `80` so we can access the Rtudio server remotely. Note you need to store a password in a variable; it will be required below for the web login.

You just need to open a web browser and point it to `localhost` if you are running Docker on your machine, or `<Your VM's IP Address>` if you are running on a cloud service.

You should see a prompt for credentials, with user defaulting to `rstudio`, and password..

Now you can run the same analysis from the RStudio console:

```
> system("Rscript readings-06.R --mean data/inflammation-*.csv")
```
{: .r}

Once you're done, stop the container with:

```
$ docker stop rstudio
```
{: .bash}


### Running a scripted R workflow on HPC with Shifter ###

We can run the same analysis on HPC through command line using Shifter. 

To get started let's pull the required R container image:

```
$ module load shifter
$ sg $PAWSEY_PROJECT -c 'shifter pull rocker/tidyverse:3.5'
```
{: .bash}

Now let's change directory to either `$MYSCRATCH` or `$MYGROUP`, e.g.

```
$ cd $MYSCRATCH
```
{: .bash}

Let's create a dedicated directory and download the sample data:

```
$ mkdir r_example
$ cd r_example
$ wget http://swcarpentry.github.io/r-novice-inflammation/data/r-novice-inflammation-data.zip
$ unzip -q r-novice-inflammation-data.zip
```
{: .bash}

With your favourite text editor, create the R file `readings-06.R` (see contents above),

and then create a SLURM script, we'll call it `rscript.sh` (remember to specify your Pawsey project ID in the script!):

```
#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --export=NONE
#SBATCH --job-name=rstudio

module load shifter

# run R script
srun --export=all shifter run rocker/tidyverse:3.5 Rscript readings-06.R --mean data/inflammation-*.csv
```
{: .bash}

Let's submit the script via SLURM:

```
$ sbatch --reservation <your-pawsey-reservation> rscript.sh
```
{: .bash}
