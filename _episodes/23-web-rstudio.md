---
title: "GUI enabled web applications: RStudio in a container"
teaching: 10
exercises: 10
questions:
objectives:
- Run interactive GUI web sessions from a container
- Setup a long running web service from a container
keypoints:
- An interactive web session can be executed as any other containerised application, via `singularity exec`
- Use the `%startscript` section of a def file to configure an image for long running services
- Launch/shutdown long running services in the background with `singularity instance start/stop`
---


### R and RStudio images

R is a popular language in several domains of science, mostly because of its statistical packages.  In particular it is nowadays highly common in data science and bioinformatics.

The group [Rocker](https://hub.docker.com/r/rocker) has published a set of R images we can use, including an Rstudio image.  To begin, let's cd into the appropriate directory:

```
$ cd $TUTO/demos/rstudio
```
{: .bash}


> ## Pull the container
>
> We want to use a [Tidyverse](https://www.tidyverse.org) container image (contains R, RStudio, data science packages).  Can you pull the  `rocker/tidyverse:3.6.1` from Docker Hub?
>
> > ## Solution
> >
> > ```
> > $ singularity pull docker://rocker/tidyverse:3.6.1
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


### Running a scripted R workflow on the shell

To begin with, we are going to run a minimalistic example taken from the workshop [Programming with R](http://swcarpentry.github.io/r-novice-inflammation/) by the Software Carpentry.  In particular, their [Episode 5](http://swcarpentry.github.io/r-novice-inflammation/05-cmdline/index.html) is the source for the dataset and the R script file; the latter has been adapted for this workshop.

Let us start with running the R script through the R container; we're going to compute average values in this example:

```
$ Rscript readings-density.R --mean inflammation-density.png data/inflammation-*.csv
```
{: .bash}


> ## Run the R script with Singularity
> 
> How would you run the above command using Singularity and the R image you just puled?
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec tidyverse_3.6.1.sif Rscript readings-density.R --mean inflammation-density.png data/inflammation-*.csv
> > ```
> > {: .bash}
> > 
> > ```
> > 5.45
> > 5.425
> > 6.1
> > 
> > [..]
> > 
> > 6.875
> > 6.025
> > 6.9
> > Saving 7 x 7 in image
> > ```
> > {: .output}
> > 
> > We even got a nice plot file out of the analysis, `inflammation-density.png`. 
> {: .solution}
{: .challenge}

So, what if want to run our R workflow using the RStudio web-based GUI interface?


### Run an interactive RStudio session

Beside developing R container images, Rocker has also some useful documentation on [Running Rocker R container with Singularity](https://www.rocker-project.org/use/singularity/).  Let's set this up together.

The Rocker documentation page suggests that an appropriate command to spawn the RStudio server is (do not run it, yet, we'll do it from the container):

```
$ rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=0  --auth-pam-helper-path=pam-helper
```
{: .bash}

Here, we're saying we want the web server to listen to port 8787 on any IP address (0.0.0.0), and then we're using another two flags to configure the authenticator.


> ## Communication ports
>
> In order to be able to use the web server, you need to ensure that the machine you are running Singularity from has opened the communication port you're using, in this case `8787`.  
> On cloud platforms, such as Nimbus at Pawsey, this will typically involve some setup in the system dashboard.  
{: .callout}


Do we need more? Yes, we need to ensure we know the username and password for authenticating; *rserver* will configure them based on the values of the container environment variables `USER` and `PASSWORD`; normally we would pick a random string for the latter.  Today we'll use "rstudiopassword".

```
$ export SINGULARITYENV_USER=$USER
$ export SINGULARITYENV_PASSWORD=rstudiopassword
$ echo $SINGULARITYENV_USER && echo $SINGULARITYENV_PASSWORD
```
{: .bash}

Lastly, containers are read-only, but RStudio will want to be able to write configuration and temporary files in the home.  Let us bind mount the current work directory as the container home (we'll use the `-B` flag).  
There's a little caveat here, due to the way the **Rocker** developers designed the container image in the recipe file.  Depending on your user ID in the host machine, the actual username in the RStudio server might be `rstudio` (*ID = 1000*) or the same as your user (otherwise); this has an impact on the container home path.  
A simple solution to this is to bind mount the fake home in the host twice, once mapped to `/home/rstudio` and once to `$HOME`. Singularity let us do this, and will handle the two bind mounts correctly.

<!-- This solution is sooo elegant, but also an over-complication to the audience -->
<!--
Let us code these conditions as follows, using a bit of bash syntax:
```
$ export HOME_USER=$USER && [ "$(id -u)" == "1000" ] && export HOME_USER=rstudio
```
{: .bash}
-->

Now we have everything we need to put together the Singularity idiomatic way to launch an interactive RStudio web server:

```
$ export SINGULARITYENV_USER=$USER
$ export SINGULARITYENV_PASSWORD=rstudiopassword
$ echo $SINGULARITYENV_USER && echo $SINGULARITYENV_PASSWORD

$ singularity exec \
    -C \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    tidyverse_3.6.1.sif \
    rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=0 --auth-pam-helper-path=pam-helper
```
{: .bash}

Note the `-C` flag for `singularity exec`, used to isolate the container from the host, including the use of a volatile `/tmp` directory instead of the host one, to better clean up the session upon exit.  As a by product, shell environment is also isolated, which is why we're defining `USER` and `PASSWORD` prefixing them with `SINGULARITYENV_`.  
If everything is fine, no output will be printed.

Now, open your web browser, and type the following as URL: `<Singularity machine IP Address>:8787`.  The `IP` can be replaced with `localhost` if you're running locally on your laptop or workstation.  
To fill the credential fields, Use `$USER` and `$PASSWORD` as printed by the commands above.  
Then we can use RStudio!

In the R console, submit the analysis script we ran earlier on from the shell:

```
> source("readings-density.R")
```
{: .r}

If you have a look at the bottom right panel, you can see some outputs files are generated, including `interactive.png`.  Click on it, and you'll get to visualise the resulting plot!

Once you're done, click on the power icon on the top right to close the session, then go back to the shell and kill the container with `Ctrl-C`.

As a final remark, note that the setup we just described could be adapted for use from a compute node in a HPC system, too, by using the HPC scheduler.


### Setup a long running RStudio web server

The procedure we just described can be convenient for interactive sessions of relatively short duration.  On the other hand, if we wanted to deploy a long running RStudio server, having to keep the terminal open isn't really handy.

Singularity has features to run containers in background.  To this end we're going to explore the subcommands of `singularity instance`.

If we need an image to be run as a background instance with Singularity, this needs to be built with a special section in the def file, namely `%startscript`.  Commands in this section are executed when the `instance` is started in background.  If no such section is provided, by default a shell will be executed, which is a bit useless for our RStudio server.

Building on the experience in the past paragraph, let us design a def file for the purpose (see `tidyverse_long.def` in the demo dir:

```
Bootstrap: docker
From: rocker/tidyverse:3.6.1

%labels
  Author Pawsey Supercomputing Centre
  Version 0.0.1

%startscript
  export R_PORT=${R_PORT:-"8787"}
  export R_ADDRESS=${R_ADDRESS:-"0.0.0.0"}

  rserver --www-port $R_PORT --www-address $R_ADDRESS --auth-none=0 --auth-pam-helper-path=pam-helper
```
{: .source}

Basically, we're starting from the `tidyverse` Docker image we used above, and then adding some commands under the `%startscript` header.  In particular, we're adding some flexibility to the `rserver <..>` command we used above, allowing for port and address to be redefined by the user through environment variables, and at the same time by providing sensible defaults.


> ## Build an image to run a RStudio instance
>
> How would you build an image called `tidyverse_long.sif`, starting from this def file?
>
> > ## Solution
> >
> > ```
> > $ sudo singularity build tidyverse_long.sif tidyverse_long.def
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

Once the container image is built, let's use it to start an instance via `singularity instance start`.  Note how the other options are the same as for the interactive session above; the only addition is the specification of a name for the instance, `myserver` in this case, that has to follow the image name:

```
$ export SINGULARITYENV_USER=$USER
$ export SINGULARITYENV_PASSWORD=rstudiopassword
$ echo $SINGULARITYENV_USER && echo $SINGULARITYENV_PASSWORD

$ singularity instance start \
    -C \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    tidyverse_long.sif \
    myserver
```
{: .bash}

```
INFO:    instance started successfully
```
{: .output}

We can check on the running instances with

```
$ singularity instance list
```
{: .bash}

```
INSTANCE NAME    PID      IMAGE
myserver         18080    /home/ubuntu/singularity-containers/demos/rstudio/tidyverse_long.sif
```
{: .output}

Note that we can run commands from the instance by referring to it as `instance://<INSTANCE-NAME>`, *e.g.*

```
$ singularity exec instance://myserver bash -c 'echo $USER $PASSWORD'
```
{: .bash}

Once we've finished with RStudio, we can shutdown the instance with

```
$ singularity instance stop myserver
```
{: .bash}

```
Stopping myserver instance of /home/ubuntu/singularity-containers/demos/08_rstudio/tidyverse_long.sif (PID=18080)
```
{: .output}
