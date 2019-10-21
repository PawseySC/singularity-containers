---
title: "Long running services with Docker"
teaching: 10
exercises: 5
questions:
objectives:
- Learn how to start containers for a running (web) service
keypoints:
- "You-ve learned how to run long-running services (like a web server) through containers"
- "Use the flag `-d` to run the containers in background"
- "Use the flag `-p <host port:<container port>` to map communication ports"
- "Additional options to manage and query containers include `--name` and `docker logs`"
---

### Starting a long-running service in a container ###
Containers are useful for running services, like web-servers etc. Many come packaged from the developers, so you can start one easily, but first you need to find the one you want to run. You can either search on [Docker Hub](https://hub.docker.com), or you can use the `docker search` command. Nginx is a popular web-server, let's look for that:

```
$ docker search nginx
```
{: .bash}

```
NAME                      DESCRIPTION                                     STARS     OFFICIAL   AUTOMATED
nginx                     Official build of Nginx.                        4719      [OK]       
jwilder/nginx-proxy       Automated Nginx reverse proxy for docker c...   877                  [OK]
richarvey/nginx-php-fpm   Container running Nginx + PHP-FPM capable ...   311                  [OK]
million12/nginx-php       Nginx + PHP-FPM 5.5, 5.6, 7.0 (NG), CentOS...   76                   [OK]
webdevops/php-nginx       Nginx with PHP-FPM                              63                   [OK]
maxexcloo/nginx-php       Framework container with nginx and PHP-FPM...   58                   [OK]
bitnami/nginx             Bitnami nginx Docker Image                      20                   [OK]
gists/nginx               Nginx on Alpine                                 8                    [OK]
evild/alpine-nginx        Minimalistic Docker image with Nginx            8                    [OK]
million12/nginx           Nginx: extensible, nicely tuned for better...   8                    [OK]
maxexcloo/nginx           Framework container with nginx installed.       7                    [OK]
webdevops/nginx           Nginx container                                 7                    [OK]
1science/nginx            Nginx Docker images based on Alpine Linux       4                    [OK]
ixbox/nginx               Nginx on Alpine Linux.                          3                    [OK]
drupaldocker/nginx        NGINX for Drupal                                3                    [OK]
yfix/nginx                Yfix own build of the nginx-extras package      2                    [OK]
frekele/nginx             docker run --rm --name nginx -p 80:80 -p 4...   2                    [OK]
servivum/nginx            Nginx Docker Image with Useful Tools            2                    [OK]
dock0/nginx               Arch container running nginx                    2                    [OK]
blacklabelops/nginx       Dockerized Nginx Reverse Proxy Server.          2                    [OK]
xataz/nginx               Light nginx image                               2                    [OK]
radial/nginx              Spoke container for Nginx, a high performa...   1                    [OK]
tozd/nginx                Dockerized nginx.                               1                    [OK]
c4tech/nginx              Several nginx images for web applications.      0                    [OK]
unblibraries/nginx        Baseline non-PHP nginx container                0                    [OK]
```
{: .output}

The official build of Nginx seems to be very popular, let's go with that:

```
$ docker pull nginx
```
{: .bash}

```
Using default tag: latest
latest: Pulling from library/nginx
386a066cd84a: Pull complete
386dc9762af9: Pull complete
d685e39ac8a4: Pull complete
Digest: sha256:3861a20a81e4ba699859fe0724dc6afb2ce82d21cd1ddc27fff6ec76e4c2824e
Status: Downloaded newer image for nginx:latest
```
{: .output}

Now, in order to run a web server such as a Nginx, we are going to use some additional Docker options, to:
* open up communication ports
* run the container in background
* give the container a specific name

We will also use a number of new Docker commands. Let's start with a known one:

```
$ docker run -p 80:80 --name=nginx nginx
```
{: .bash}

The option `-p 80:80` tells Docker to map port `80` on the host to port `80` in the container, so you can communicate with it.

We've also introduced a second new option: `--name`. Docker automatically names containers, but these names don't reflect the purpose of the container (e.g. `pensive_booth` was the name of the nginx container I ran for this example without the `--name` option).  You can name your container whatever you want, but it's helpful to give it a name similar to the container you're using, or the specific service or application you're running in the container. This can be useful when we need to act on the container while it's running, e.g. to stop it, or to get logs, as wwe'll see soon. In this example, we've called our container `nginx`. 

Note also that we didn't tell docker what program to run, that's baked into the container in this case. More on that later.

Now, go to your browser and in the address bar enter `localhost` if you are running Docker on your machine, or `<Your VM's IP Address>` if you are running on a cloud service. You should see a page with a **Welcome to nginx!** message. On your terminal where you ran the docker command, you'll see some log information:

```
172.17.0.1 - - [30/Nov/2016:18:07:59 +0000] "GET / HTTP/1.1" 200 612 "-" "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.11; rv:49.0) Gecko/20100101 Firefox/49.0" "-"
2016/11/30 18:07:59 [error] 7#7: *1 open() "/usr/share/nginx/html/favicon.ico" failed (2: No such file or directory), client: 172.17.0.1, server: localhost, request: "GET /favicon.ico HTTP/1.1", host: "localhost"
172.17.0.1 - - [30/Nov/2016:18:07:59 +0000] "GET /favicon.ico HTTP/1.1" 404 169 "-" "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.11; rv:49.0) Gecko/20100101 Firefox/49.0" "-"
172.17.0.1 - - [30/Nov/2016:18:07:59 +0000] "GET /favicon.ico HTTP/1.1" 404 169 "-" "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.11; rv:49.0) Gecko/20100101 Firefox/49.0" "-"
2016/11/30 18:07:59 [error] 7#7: *1 open() "/usr/share/nginx/html/favicon.ico" failed (2: No such file or directory), client: 172.17.0.1, server: localhost, request: "GET /favicon.ico HTTP/1.1", host: "localhost"
```
{: .output}

That's a good start, but you now have a terminal tied up with nginx, and if you hit `CTRL-C` in that terminal, your web-server dies. 

To practice a different solution, let's first stop this container. How to? Use the `docker stop` command!

```
$ docker stop nginx
```
{: .bash}

```
nginx
```
{: .output}

Check that it's gone:

```
$ docker ps
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
```
{: .output}

Now let's go ahead with using the Docker option `-d` to run the container in the background instead (daemon mode):

```
$ docker run -d -p 80:80 --name=nginx nginx
```
{: .bash}

```
48a2dca14407484ca4e7f564d6e8c226d8fdd8441e5196577b2942383b251106
```
{: .output}

Go back to your browser, reload `localhost` (or `<Your VM's IP Address>`), and you should get the page loaded again.

We can view the logs of our nginx service with the `docker logs` command, followed by the container name, `nginx`:

```
$ docker logs --follow nginx
```
{: .bash}

```
172.17.0.1 - - [30/Nov/2016:18:18:40 +0000] "GET / HTTP/1.1" 304 0 "-" "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.11; rv:49.0) Gecko/20100101 Firefox/49.0" "-"
```
{: .output}

This gives us a live look at what is going on in our nginx container (try reloading your webpage and see what the logs look like).  Note that the `--follow` option keeps the terminal open and will update the logs in real-time.  If you omit it, `docker logs` would simply display the last few lines of the log file.

If you hit `CTRL-C` now, your container is still running, in the background. You can see this with the `docker ps` command:

```
$ docker ps
```
{: .bash}

```
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                           NAMES
48a2dca14407        nginx               "nginx -g 'daemon off"   3 minutes ago       Up 3 minutes        443/tcp, 0.0.0.0:80->80/tcp   nginx
```
{: .output}

You can open a shell into the running container, if you wish, using `docker exec`. This can be handy for debugging:

```
$ docker exec -t -i nginx /bin/bash
```
{: .bash}

```
root@48a2dca14407:/# 
```
{: .output}

```
root@48a2dca14407:/# exit   # or hit CTRL-D
```
{: .bash}

<!-- does not work with latest nginx images.. ps is not installed!
Let's take a look at the processes running in our nginx server:
```
$ docker exec -t -i nginx /bin/bash
```
{: .bash}

```
root@48a2dca14407:/# ps auxww | grep ngin                                                                                                                
```
{: .bash}

```
root         1  0.0  0.0  31764  5164 ?        Ss   18:58   0:00 nginx: master process nginx -g daemon off;
nginx        7  0.0  0.0  32156  2900 ?        S    18:58   0:00 nginx: worker process
root        15  0.0  0.0  11128  1036 ?        S+   18:59   0:00 grep ngin
```
{: .output}
-->

Finally, let's stop this container:

```
$ docker stop nginx
```
{: .bash}

```
nginx
```
{: .output}


> ## Start a Jupyter web server using a container ##
> 
> How would you use the container image `jupyter/scipy-notebook` to start a Jupyter webserver?
> 
> A couple of extra requirements:
> 
> * assign the name `jupyter` to the container
> * map the current host working directory as `/home/jovyan/data` (see previous episode)
> 
> Some hints:
> 
> * by default, that image wants to use port `8888` for the container
> * use the default HTTP port `80` for the host, to navigate to the corresponding page in the web browser using default syntax
> * no command needs to be specified for that container, the default behaviour is to start the webserver
> 
> After you have started it, use `docker logs` to look for an access `token` string.
> 
> Once you are done, stop the container.
> 
> > ## Solution ##
> > 
> > Pull the container:
> > 
> > ```
> > $ docker pull jupyter/scipy-notebook
> > ```
> > {: .bash}
> > 
> > Start the webserver:
> > 
> > ```
> > $ docker run -d -p 80:8888 --name=jupyter -v `pwd`:/home/jovyan/data jupyter/scipy-notebook
> > ```
> > {: .bash}
> > 
> > Inspect the logs:
> > 
> > ```
> > $ docker logs jupyter
> > ```
> > {: .bash}
> > 
> > Use your web browser to go to `localhost` if you are running Docker on your machine, or `<Your VM's IP Address>` if you are running on a cloud service.
> > 
> > Stop the container:
> > 
> > ```
> > $ docker stop jupyter
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}
