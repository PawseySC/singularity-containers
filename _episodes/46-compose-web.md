---
title: "Run multi-component web services with Docker Compose"
teaching: 5
exercises: 10
questions:
objectives:
- Get an idea of how Docker Compose allows to setup multi-component web services
keypoints:
- Docker Compose uses a YAML file to layout container specifications, including interactions among containers
- Use `docker-compose up -d` to start the service, and `docker-compose down` to shut it down
- Inspect the shell output of the running service via `docker-compose logs`
---


### Why Docker Compose?

In previous episodes we've seen how to use Singularity to spawn interactive web sessions, or even long running web services.  Why would we use another tool for this?

Well, in this context [Docker Compose](https://docs.docker.com/compose/) provides a key additional capability compared to Singularity.  
Building on top of Docker, it allows for setting up, launching and coordinating multiple containers at once.  This is extremely powerful for any web service that requires multiple components in order to work properly.

Typical examples of multi-component services, that can be of relevance to research communities, are:
* services requiring a dedicated database, such as lab notebooks, or collaboration platforms (*i.e.* the example we'll run here)
* any service that you need to secure through a reverse proxy and the HTTPS protocol (literally any, including RStudio, Jupyter Hub...)

Making use of Docker, Compose cannot be run on shared HPC systems right now, only on a local laptop/workstation and on Cloud.


### A sample specification file

The example here is taken from the [deployment guide](https://hackmd.io/c/codimd-documentation/%2Fs%2Fcodimd-docker-deployment) of the open source collaboration platform [CodiMD](https://hackmd.io/c/codimd-documentation).

Specification files for Docker Compose are in YAML format and are normally called `docker-compose.yml`:

```
version: "3"
services:
  database:
    image: postgres:11.6-alpine
    environment:
      - POSTGRES_USER=codimd
      - POSTGRES_PASSWORD=Chang3Th1sPassw0rd
      - POSTGRES_DB=codimd
    volumes:
      - database-data:/var/lib/postgresql/data
    restart: always
  codimd:
    image: nabo.codimd.dev/hackmdio/hackmd:2.1.0-cjk
    environment:
      - CMD_DB_URL=postgres://codimd:Chang3Th1sPassw0rd@database/codimd
      - CMD_USECDN=false
    depends_on:
      - database
    ports:
      - 3000:3000
    volumes:
      - upload-data:/data/codimd_test/hackmd/app/public/uploads
    restart: always
volumes:
  database-data: {}
  upload-data: {}
```
{: .source}

We don't want to understand the entire file, just point out some aspects:

* the file specifies *two* service containers, one called `database` and one called `codimd`;

* for each container, typical container properties are defined, including image name, attached volumes, communication ports, environment variables;

* note how dependencies between containers can be specified with the keyword `depends_on`.

Before moving on, let us also notice that the `codimd` container exposes binds to the communication port `3000` in the host.


### Run a web service with Docker Compose

Let's cd into the appropriate demo directory, where a copy of the YAML file above is:

```
$ cd $TUTO/demos/codimd
$ ls 
```
{: .bash}

```
docker-compose.yml
```
{: .output}

We can launch the web service described in this YAML file by means of (`-d` is for *daemon* mode, *i.e.* it will run the processes in background):

```
$ docker-compose up -d
```
{: .bash}

> ## Communication ports
>
> In order to be able to use the web server, we need to ensure that the machine we're running Docker Compose from has opened the communication port we need, in this case `3000`.  
> On cloud platforms, such as Nimbus at Pawsey, this will typically involve some setup in the system dashboard.  
{: .callout}

Now, let's open our web browser, and type the following as URL: `<Docker Compose machine IP Address>:3000`.  The `IP` can be replaced with `localhost` when running locally on a laptop or workstation.  

We'll get to the home page of CodiMD, where we can create a user, login, and then start creating notes in Markdown, ready to be shared with the world!

To finish off, there are a couple more of handy commands to mention.

If we need to have a look at the logs of the web services, we can use the command:

```
$ docker-compose logs -f
```
{: .bash}

and then hit `Ctrl-C` when we're done watching.

To shut down our web service we can use:

```
$ docker-compose down
```
{: .bash}

This command will still preserve the volumes containing the data related to the services (users, notes,...).  If we need to get rid of the volumes, too, all we need to do is to add the flag `-v`.

To resume the service, just execute `docker-compose up` again.
