---
title: "Making Python not awful with containers"
teaching: 0
exercises: 20
questions:
objectives:
- Run Python applications using container
- Use JupyterHub in a container
keypoints:
- Containers are great way to manage Python workflows
---

### Why can Python be Awful? ###

Python is a great language for doing all kinds of work, but sometimes it can present issues...

![Python can get messy]({{ page.root }}/fig/python-complexity-cartoon.png)


### Build and run a basic Python app ###

Let's start by running a very simple hello world example with a basic Python container.  To start, `cd` 
to the `/demos/07_python` directory and create a basic `helloworld.py`:

```
#!/usr/bin/env python

print("Hello World")
```
{: .source}

And we can run this without even building a new Docker image:

```
chmod +x helloworld.py

docker run --rm -v $(pwd):/app -w /app python:3.6-slim python ./helloworld.py
```
{: .bash}


```
Unable to find image 'python:3.6-slim' locally
3.6-slim: Pulling from library/python
fc7181108d40: Pull complete
01f3b6f74410: Pull complete
d8e6069e053d: Pull complete
e6bd64ad4fe3: Pull complete
a034e22b86de: Pull complete
Digest: sha256:0d095570901e7cf0dac93ba4d0ee75ee2364715338b0396874589c9cc23343b6
Status: Downloaded newer image for python:3.6-slim

Hello World
```
{: .output}


### Adding Python modules to a container ###

We want to be able leverate other Python modules to do some acutal work, so we need to be able
to build Python containers that use tools like `pip` and `conda`.

For this example, `cd` to the `demos/07_python/logistic-regression` and we'll trying installing some
ML packages to run a logistic regression model.

We have both a Jupyter notebook (`LogisticRegression.ipynb`) and source code (`logreg.py`) that we'll use.
An example Dockerfile is also provided.

> Build a custom Jupyter/Python container
> For this example you should try building your own Docker container.  Some hints:
>
> * Decide if you should use a base image (the [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html) might be a good place to start)
> * You'll need to install the Plotly pacakge
>
> > ## Solution ##
> >
> > ```
> > FROM jupyter/datascience-notebook:latest
> >
> > RUN conda install plotly
> > ```
> > {: .source}
> >
> > You can build it with
> >
> > ```
> > docker build -t logreg-py .
> > ```
> > {: .bash}
> {: .solution}

After you've built your image we'll start up the container and log into our Jupyter notebook server:

```
docker run --name=jupyter-logreg -d -p 80:8888 -v $(pwd):/home/jovyan -w /home/jovyan jupyter-logreg:latest
```
{: .bash}

Then you'll need to use `docker logs` to find out the access key for your Jupyter server:

```
docker logs jupyter-logreg
```
{: .bash}

```
Executing the command: jupyter notebook
[I 02:55:14.895 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
[I 02:55:16.361 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.7/site-packages/jupyterlab
[I 02:55:16.361 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
[I 02:55:16.363 NotebookApp] Serving notebooks from local directory: /home/jovyan
[I 02:55:16.363 NotebookApp] The Jupyter Notebook is running at:
[I 02:55:16.363 NotebookApp] http://(54bdc6f5c9cd or 127.0.0.1):8888/?token=82be2f08380f138a92ad9c6323f7fcfd124c453144492d3f
[I 02:55:16.363 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 02:55:16.387 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/jovyan/.local/share/jupyter/runtime/nbserver-6-open.html
    Or copy and paste one of these URLs:
        http://(54bdc6f5c9cd or 127.0.0.1)/?token=82be2f08380f138a92ad9c6323f7fcfd124c453144492d3f
```
{: .output}

Point your web browser to **http://<nimbus-ip/?token=<token>** and we'll go through the notebook.

### Running Python code without Jupyter ###

Because our Jupyter container has a functioning version of Python installed, we can also run our container
and Python code directly from the command line:

```
docker run --name=python-logreg -v $(pwd):/home/jovyan -w /home/jovyan jupyter-logreg:latest python ./logreg.py
```

You should have plot now saved as a `.png` file.

### Machine Learning Example ###

Let's try another ML example.  We'll build another Docker image that uses Jupyter.  For this  this example, `cd` to the `demos/07_python/image-classification`.

In here we have some image data (`dataset/`), a Jupyter notebook, and a Dockerfile.

Just like before, we'll have you build your own Jupyter container.

> Build a custom Jupyter/Python container
> For this example you should try building your own Docker container.  Some hints:
>
> * Decide if you should use a base image (the [Jupyter Docker Stacks](https://jupyter-docke
stacks.readthedocs.io/en/latest/using/selecting.html) might be a good place to start)
> * You'll need to install 2 modules: **mahotas** annd **opencv-python**
>
> > ## Solution ##
> >
> > ```
> > FROM jupyter/datascience-notebook:latest
> >
> > RUN conda install mahotas
> > RUN pip install opencv-python
> > ```
> > {: .source}
> >
> > You can build it with
> >
> > ```
> > docker build -t jupyter-image .
> > ```
> > {: .bash}
> {: .solution}

After you've built your image we'll start up the container and log into our Jupyter notebook
server:

```
docker run --name=jupyter-image -d -p 80:8888 -v $(pwd):/home/jovyan -w /home/jovyan jupyter-image:latest
```
{: .bash}

Then you'll need to use `docker logs` to find out the access key for your Jupyter server:

```
docker logs jupyter-images
```
{: .bash}

```
Executing the command: jupyter notebook
[I 02:55:14.895 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
[I 02:55:16.361 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.7/site-packages/jupyterlab
[I 02:55:16.361 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
[I 02:55:16.363 NotebookApp] Serving notebooks from local directory: /home/jovyan
[I 02:55:16.363 NotebookApp] The Jupyter Notebook is running at:
[I 02:55:16.363 NotebookApp] http://(54bdc6f5c9cd or 127.0.0.1):8888/?token=82be2f08380f138a92ad9c6323f7fcfd124c453144492d3f
[I 02:55:16.363 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 02:55:16.387 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/jovyan/.local/share/jupyter/runtime/nbserver-6-open.html
    Or copy and paste one of these URLs:
        http://(54bdc6f5c9cd or 127.0.0.1)/?token=82be2f08380f138a92ad9c6323f7fcfd124c453144492d3f
```
{: .output}

Point your web browser to **http://<nimbus-ip/?token=<token>** and we'll go through the notebook.
