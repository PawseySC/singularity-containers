---
title: "GUI applications using X11 windows: the Gnuplot example"
teaching: 10
exercises: 5
questions:
objectives:
- Run an X11 enabled GUI application from a container
keypoints:
- If launching an X11 application remotely, you need to enable X11 forwarding for the SSH connection
- Build the image for an X11 application including the package `xauth`, for compatibility with the Docker runtime
- When running the X11 application, bind mount the `Xauthority` file
---


### Graphical applications using the X Window system

A good number of Linux scientific packages provide a graphical interface by means of desktop windows.  The *X Window System*, or *X11*, is one of the most common of this kind.  It allows an application to open graphical windows not only locally, but also remotely using a mechanism called *X11 forwarding.  In addition, dedicated clients exist for macOS and Windows that enable remote forwarding on this platforms.

In this episode we're going to learn how to run a X11 enabled application from inside a container.  As an example, we're using Gnuplot, a simple yet powerful plotting utility.


### Requirements for the host system (local or remote)

Regardless whether we're need to run our X11 application locally or remotely, our local machine needs an X11 system.  Linux boxes should all be good to go.  On macOS, you'll need to install [XQuartz](https://www.xquartz.org).  On Windows, there are a few alternatives; [Cygwin/X](https://x.cygwin.com) seems quite common.  
If we're running our application by remotely connecting to a Linux machine, that machine will need X11 setup, too.  Moreover, we'll need to enable X11 forwarding when establishing the SSH connection to that machine.  This is done via either

```
$ ssh -X
```
{: .bash}

or 

```
$ ssh -Y
```
{: .bash}

depending on whether we're connecting from Linux or macOS; in Windows, the graphical SSH client will have its own dedicated setting for X11 forwarding, that needs to be enabled.


### Extra package in the container image

There is one small requirement for the image that packages the X11 enabled application.  It's a utility called `xauth`, used to authenticate the X client in the container towards the X server of the local machine, so that the latter can go ahead in opening windows and receiving commands from the former.  
In practice, the Dockerfile (or def file) for the image needs to contain an `apt` command to install `xauth`.  *E.g.* in a Dockerfile:

```
[..]

RUN apt-get install -y xauth

[..]
```
{: .source}

To be precise, this additional package is not required when running X11 containerised applications with Singularity.  However, it would be required when running it with Docker, so it's good practice to embed it in order to enforce compatibility across different container engines.


### Extra flag when running the container with Singularity

One last thing to know to run our example is that we need to bind mount an *Xauthority* file in the container.  This file contains a secret string, or *magic cookie*, that is used by the X client to authenticate with the X server.  

Usually this file is located in the user's home directory, so the mounting will look like:

```
-B ~/.Xauthority
```
{: .bash}

Less often, a dedicated location is used for this file which is contained in the `$XAUTHORITY` variable.  In this case, this will do the job: `-B $XAUTHORITY`.


### Let's run a Gnuplot script and *see* what happens!

First remember that, in order for this example to work as intended, your local machine needs to have a X server installed (see above).  
Also, if you are logging in to a remote machine, you need to enable X11 forwarding using `ssh -X` (from Linux), or `ssh -Y` (from macOS), or from Windows whatever SSH client setting is required.  

Let's `cd` into the appropriate directory:

```
$ cd $TUTO/demos/gnuplot
```
{: .bash}

The sample script we're going to use comes from the Gnuplot Collection of Demo Scripts (see bottom of [this page](http://gnuplot.sourceforge.net/demo/hidden2.html)).  
We're going to use the image `docker://marcodelapierre/gnuplot:5.2.2_4`.  Let's start `gnuplot` from this image.  Don't forget the additional bind mounting:

```
$ singularity exec -B ~/.Xauthority docker://marcodelapierre/gnuplot:5.2.2_4 gnuplot
```
{: .bash}

```
gnuplot>
```
{: .output}

Now, from the gnuplot shell, execute the following:

```
gnuplot> load('tori.gp')
```
{: .source}

Whooa, we got doughnuts!

When you're done, close the X11 window by hitting `q`, and then close the gnuplot shell using `exit`.
