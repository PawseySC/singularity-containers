---
title: "Build and share your own container image"
teaching: 10
exercises: 10
questions:
objectives:
- Build a container image
- Share an image with others
keypoints:
- Build images using `sudo singularity build`
- Use the remote builder with the flag `-r`, if you need to build images from a machine where you don't have sudo rights
- You can share your Singularity Image File with others, as you would do with any other (big) file
- Upload images to a web registry with `singularity push` (Sylabs account required)
---


### Building a basic container

Singularity can build container images in different formats.  Let's focus on the Singularity Image Format, *SIF*, which is the one typically adopted to ship production-ready containers.  
This example is adapted from this well crafted [Singularity Tutorial](https://github.com/ArangoGutierrez/Singularity-tutorial).

Let us cd into the appropriate directory:

```
$ cd $TUTO/demos/lolcow
```
{: .bash}

To build an image we need a recipe, or **definition file**, in the Singularity language.  You'll learn more on how to write one in a dedicated episode.

Here, let's use one, `lolcow.def`, to build our first image.  To this end we're using `sudo singularity build`, followed by the filename (and path) we decide to attribute to the container image, and then by the filename of the def file to be used:

```
$ sudo singularity build lolcow.sif lolcow.def
```
{: .bash}

This build should take approximately 5 minutes to complete:

```
INFO:    Starting build...
[..]
INFO:    Running post scriptlet
[..]
INFO:    Adding help info
INFO:    Adding labels
INFO:    Adding environment to container
INFO:    Adding runscript
INFO:    Creating SIF file...
INFO:    Build complete: lolcow.sif
```
{: .output}

At the end, you'll find the final image file, `lolcow.sif`, in your directory:

```
$ ls
```
{: .bash}

```
lolcow.def lolcow.sif
```
{: .output}


### *Sudo* privileges with Singularity

As the image builds, let's discuss this important matter.  

Singularity does not allow for privileges escalation.  
In other words, if you are a standard user and you run `singularity`, any command inside the container will be run with the privileges of the standard user, *i.e.* without admin powers.  If you try and `sudo` from inside the container you will get an error.  
On the other hand, if your user can run with *sudo*, and if you then decide to run Singularity as `sudo singularity`, then you will run any command from inside the container with admin powers.  
This design contributes to making Singularity safe to run on HPC: users without admin rights are unable to escalate their privileges from inside the containers.

However, when building a container image you might need to install software using commands that require admin rights, *e.g.* `apt get` in Ubuntu/Debian or `yum` in Centos.  To achieve this, you need to run `sudo singularity build`, as we have just done above.  
This requirement implies that you must carry out your build in a machine where you DO have admin rights. Ruling out HPC systems, suitable platforms to build container images can be your laptop, an office workstation, or a virtual machine on the cloud.


### Remote build

What if you need to build an image from a system where you don't have admin privileges, *i.e.* you can't run commands with *sudo*?

Singularity offers the option to run a build remotely, using the **Sylabs Remote Builder**; once again (see below) you will need a Sylabs account and an `access token` to use this feature.  If this is the case, just use `singularity build -r` to proceed with the remote build.  Once finished, the image will be downloaded so that it's ready to use:

```
$ singularity build -r lolcow_remote.sif lolcow.def
```
{: .bash}

```
INFO:    Remote "default" added.
INFO:    Authenticating with remote: default
INFO:    API Key Verified!
INFO:    Remote "default" now in use.
INFO:    Starting build...
[..]
INFO:    Running post scriptlet
[..]
INFO:    Adding help info
INFO:    Adding labels
INFO:    Adding environment to container
INFO:    Adding runscript
INFO:    Creating SIF file...
INFO:    Build complete: /tmp/image-699539270
WARNING: Skipping container verifying
 67.07 MiB / 67.07 MiB  100.00% 14.18 MiB/s 4s
```
{: .output}

At the time of writing, when using the Remote Builder you won't be able to use the `%files` header in the def file, to copy host files into the image.


> ## Use the newly created container
>
> Once the image has finished building, how would you run the command `fortune` from inside this container?
>
> > ## Solution
> >
> > ```
> > $ singularity exec lolcow.sif fortune
> > ```
> > {: .bash}
> >
> > ```
> > Whenever you find that you are on the side of the majority, it is time
> > to reform.
> > 		-- Mark Twain
> > ```
> > {: .output}
> {: .solution}
> 
> Now, try and run this pipe of commands: `bash -c 'fortune | cowsay | lolcat'`.
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec lolcow.sif bash -c 'fortune | cowsay | lolcat'
> > ```
> > {: .bash}
> > 
> > You will get something similar to this, hopefully just more colourful:
> > 
> > ```
> >  _______________________________________
> > / Have a place for everything and keep  \
> > | the thing somewhere else; this is not |
> > | advice, it is merely custom.          |
> > |                                       |
> > \ -- Mark Twain                         /
> >  ---------------------------------------
> >         \   ^__^
> >          \  (oo)\_______
> >             (__)\       )\/\
> >                 ||----w |
> >                 ||     ||
> > ```
> > {: .output}
> {: .solution}
> 
> Great, we've just containerised a cow that cites novelists!  
{: .challenge}


### Share your container image

Now that you've built your container image, you might want to run it on other systems, or share it with collaborators.

The simplest way to achieve this is to remember that a SIF image is just a file, so... you can transfer it across systems using Linux command line utilities like `scp` or `rsync`, or even graphical applications such as `Filezilla`.  
Just remember that images can be quite large, typically ranging from tens of MBs up to several GBs.  The *lolcow* image we created is about 70 MB, but for instance a typical *RStudio* image is well above 1 GB.

If you want to keep the images publicly available, as they are just files you can store them in a server that is accessible through HTTP or FTP, or design a setup based upon open source image registry solutions such as [Harbor](https://goharbor.io).

Sylabs offer their own image hosting platform, the [**Sylabs Cloud Library**](https://cloud.sylabs.io), which is currently free upon signup.  Let's see how this works.  If you don't want to signup, just skip the hands-on and follow the demo.  
Once you create an account, you'll need to click on your account name on the top right of the page, select `Access Tokens`, then create a token, and copy it to the clipboard.  
Then you can configure the machine you're using for building container images, so that you can also push them to the Cloud Library:

```
$ singularity remote login
```
{: .bash}

```
Generate an API Key at https://cloud.sylabs.io/auth/tokens, and paste here:
API Key:
```
{: .output}

Now paste the token you had copied to the clipboard end press `Enter`:

```
INFO:    API Key Verified!
```
{: .output}

You are now ready to push your image to the Cloud Library, *e.g.* via `singularity push`:

```
$ export MY_SYLABS_USER="my_sylabs_user"
$ singularity push -U lolcow.sif library://$MY_SYLABS_USER/default/lolcow:30oct19
```
{: .bash}

```
WARNING: Skipping container verifying
 67.08 MiB / 67.08 MiB [==================================================================================================================================] 100.00% 6.37 MiB/s 10s
```
{: .output}

Note the use of the flag `-U` to allow pushing unsigned containers (for more information see this Singularity [documentation page](https://sylabs.io/guides/3.5/user-guide/signNverify.html)).  
Also note once again the format for the registry: `<user>/<project>/<name>:<tag>`.

Finally, you (or other peers) are now able to pull your image from the Cloud Library:

```
$ singularity pull -U library://$MY_SYLABS_USER/default/lolcow:30oct19
```
{: .bash}

```
INFO:    Downloading library image
 67.07 MiB / 67.07 MiB [===================================================================================================================================] 100.00% 8.10 MiB/s 8s
WARNING: Skipping container verification
INFO:    Download complete: lolcow_30oct19.sif
```
{: .output}
