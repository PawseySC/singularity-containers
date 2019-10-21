---
title: "Cleaning up Docker"
teaching: 8
exercises: 2
questions:
objectives:
- Learn how to remove containers and images from your machine when you no longer need them
keypoints:
- "Cleaning up containers and images is a two-step process"
- "Remove stopped containers with `docker rm`"
- "Delete unnecessary images with `docker rmi`"
- "`docker run --rm` allows to automatically remove containers at completion"
---

### Cleaning up ###

Eventually, you may want to clean out unnecessary images and the cache of containers, to reclaim space or just keep things tidy. Cleaning up involves two steps, removing the containers that you've run first, then removing the images themselves.

To remove the containers, including those that have exited and are still in the cache, use `docker ps --all` to get the IDs, then `docker rm` followed by the container ID(s) or the container name(s):

```
$ docker ps --all
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS                        PORTS               NAMES
48a2dca14407        nginx               "nginx -g 'daemon off"   10 minutes ago      Exited (0) 4 minutes ago                          pensive_booth
1bf15ca4ba89        nginx               "nginx -g 'daemon off"   20 minutes ago      Exited (0) 16 minutes ago                         amazing_hugle
97b1e86df1f1        ubuntu              "/bin/bash"              37 minutes ago      Exited (0) 37 minutes ago                         backstabbing_brown
1688f55c3418        ubuntu              "/bin/bash"              37 minutes ago      Exited (127) 20 minutes ago                       ecstatic_hugle
c69d6f8d89bd        ubuntu              "/bin/bash"              41 minutes ago      Exited (0) 41 minutes ago                         sad_keller
960588723c36        ubuntu              "/bin/echo 'hello wor"   45 minutes ago      Exited (0) 45 minutes ago                         suspicious_mestorf
ffbb0f60bda6        ubuntu              "/bin/echo 'hello wor"   51 minutes ago      Exited (0) 51 minutes ago                         mad_booth
```
{: .output}

```
$ docker rm ffbb0f60bda6 960588723c36 # cleaning by ID
```
{: .bash}

```
ffbb0f60bda6
960588723c36
```
{: .output}

```
$ docker rm sad_keller ecstatic_hugle # cleaning by name
```
{: .bash}

```
sad_keller
ecstatic_hugle
```
{: .output}

```
$ docker ps --all
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS                      PORTS               NAMES
48a2dca14407        nginx               "nginx -g 'daemon off"   23 minutes ago      Exited (0) 17 minutes ago                       pensive_booth
1bf15ca4ba89        nginx               "nginx -g 'daemon off"   33 minutes ago      Exited (0) 29 minutes ago                       amazing_hugle
97b1e86df1f1        ubuntu              "/bin/bash"              50 minutes ago      Exited (0) 50 minutes ago                       backstabbing_brown
```
{: .output}

You can construct your own one-liner to clean up everything, if you wish. Docker will refuse to delete something that's actively in use, so you won't screw things up too badly this way.

```
$ docker rm `docker ps --all -q`
```
{: .bash}

```
48a2dca14407
1bf15ca4ba89
97b1e86df1f1
```
{: .output}

Now that we've removed the containers, let's clean up the images they came from. This uses the `docker images` and `docker rmi` commands, in a similar manner:

```
$ docker images
```
{: .bash}

```
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
ubuntu              latest              4ca3a192ff2a        22 hours ago        128.2 MB
nginx               latest              abf312888d13        2 days ago          181.5 MB
```
{: .output}

```
$ docker rmi ubuntu
```
{: .bash}

```
Untagged: ubuntu:latest
Untagged: ubuntu@sha256:3b64c309deae7ab0f7dbdd42b6b326261ccd6261da5d88396439353162703fb5
Deleted: sha256:4ca3a192ff2a5b7e225e81dc006b6379c10776ed3619757a65608cb72de0a7f6
Deleted: sha256:2c2e0ef08d4988122fddadfe7e84d7b0aae246e6fa805bb6380892a325bc5216
Deleted: sha256:918fe8ae141d41d7430eb887e57a21d76fb0181317ec4a68f7abbd17caab034a
Deleted: sha256:00f434fa2fa1a0558f8740af28aef3d6ee546a8758c0a37dddee3f65b5290e7a
Deleted: sha256:f9545ee77b4b95c887dbebc6d08325be354274423112b0b66f7288a2cf7905cb
Deleted: sha256:d29d52f94ad5aa750bd76d24effaf6aeec487d530e262597763e56065a06ee67
```
{: .output}

You need to be sure to stop and remove a container before removing its image.  If not, you'll see an error about child images:

```
$ docker ps --all
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS                        PORTS                  NAMES
bf0a1726909a        nginx               "nginx -g 'daemon of…"   18 minutes ago      Exited (0) 16 minutes ago                            eloquent_raman
```
{: .output}

```
$ docker rmi nginx
```
{: .bash}

```
Error response from daemon: conflict: unable to remove repository reference "nginx" (must force) - container bf0a1726909a is using its referenced image 71c43202b8ac
```
{: .error}


> ## Clean up Miniconda containers from previous episode ##
> 
> Find and remove the stopped Miniconda containers you ran in the previous episode. Do not remove the container image, we'll use it later.
> 
> Hint: you can identify them as they will correspond to the image `continuumio/miniconda3:4.5.12`.
> 
> > ## Solution ##
> > 
> > Display stopped containers: 
> > 
> > ```
> > $ docker ps -a
> > ```
> > {: .bash}
> > 
> > Remove `miniconda3` stopped containers: identify their IDs from previous output, then 
> > 
> > ```
> > $ docker rm <list of IDs>
> > ```
> > {: .bash}
> > 
> > Do not use `docker rmi` as we don't want to remove the container image.
> {: .solution}
{: .challenge}


> ## Best practices ##
> 
> * Use `docker run` with the `--rm` flag when you know you won't want to re-start a container
> * If you use containers heavily, clean up the images from time to time
{: .callout}
