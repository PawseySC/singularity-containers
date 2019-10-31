---
title: "Making Python not awful with containers"
teaching: 15
exercises: 15
questions:
objectives:
- Run Python applications using a container
- Use Jupyter Notebook in a container
keypoints:
- "You can use idiomatic syntaxes of `singularity exec` to launch interactive sessions of Jupyter Notebooks, similar to RStudio"
- "Run the same workflow through command line or via Notebook using the same container"
---


### Why can Python be Awful?

Python is a great language for doing all kinds of work, but sometimes it can present issues...

![Python can get messy]({{ page.root }}/fig/python-complexity-cartoon.png)


### Run a basic Python script from a container

Let's start by running a very simple hello world example with a basic Python container.  To start, cd into the demo directory:

```
$ cd $SC19/demos/09_python
```
{: .bash}

There is an executable script called `helloworld.py`. 


> ## Use a Python container
> 
> Can you use Singularity to run `helloworld.py` using the container `python:3.8-slim`?
> 
> > ## Solution
> > 
> > ```
> > $ singularity exec docker://python:3.8-slim ./helloworld.py
> > ```
> > {: .bash}
> > 
> > ```
> > INFO:    Converting OCI blobs to SIF format
> > INFO:    Starting build...
> > [..]
> > INFO:    Creating SIF file...
> > INFO:    Build complete: /data/singularity/.singularity/cache/oci-tmp/c8f8cf71e80d0ca6d71858ac38dda8027637b04ade825b6f418b4e2da832c63a/python_3.8-slim.sif
> > 
> > Hello World!
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


### Build and run a Python container for machine learning

We want to be able to leverage other Python modules to do some actual work, so we need to build Python containers that use tools like `pip` and `conda`.

For this example we'll try installing some ML packages to run a logistic regression model. 


> ## Design the def file
> 
> Think of how you would write a def file that 
> * uses `jupyter/datascience-notebook:latest` from Docker Hub as base image 
> * installs the *Plotly* package via `/opt/conda/bin/conda install plotly=3.10`
> * we'll use this image as an interactive session, so do not care about `%startscript`
> 
> > ## When you are done check this solution
> > 
> > ```
> > Bootstrap: docker
> > From: jupyter/datascience-notebook:latest
> > 
> > %post
> >   /opt/conda/bin/conda install plotly=3.10
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}

Note that we need the full path to the `conda` executable when using the base image from `jupyter/`. Because of the way the Docker image was built, without the full path the binary would not be found.

Now, change directory to `logistic-regression`:

```
$ cd $SC19/demos/09_python/logistic-regression
```
{: .bash}


> ## Build the container image 
> 
> The def file above is present as `plotly.def`. Use it to build an image called `plotly.sif`. **Hint**: remember to use `sudo`!
> 
> > ## Solution
> > 
> > ```
> > $ sudo singularity build plotly.sif plotly.def
> > ```
> > {: .bash}
> {: .solution}
{: .challenge}


In the same directory we have both a Jupyter notebook (`LogisticRegression.ipynb`) and source code (`logreg.py`) that we're going to use.

After you've built your image we'll start up the container and log into our Jupyter notebook server (note the idiomatic expression for Jupyter Notebooks with Singularity):

```
$ singularity exec -C -B $(pwd):$HOME plotly.sif jupyter notebook --no-browser --port=8888 --ip 0.0.0.0 --notebook-dir=$HOME
```
{: .bash}

The output on the terminal will look like:

```
[I 14:29:20.678 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.7/site-packages/jupyterlab
[I 14:29:20.679 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
[I 14:29:20.682 NotebookApp] Serving notebooks from local directory: /home/ubuntu
[I 14:29:20.682 NotebookApp] The Jupyter Notebook is running at:
[I 14:29:20.682 NotebookApp] http://nimbus1:8888/?token=b2a37a0b893cdbb5949aaad1d27e8fcff80841b0a9ae5b1c
[I 14:29:20.682 NotebookApp]  or http://127.0.0.1:8888/?token=b2a37a0b893cdbb5949aaad1d27e8fcff80841b0a9ae5b1c
[I 14:29:20.682 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 14:29:20.688 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///home/ubuntu/.local/share/jupyter/runtime/nbserver-17-open.html
    Or copy and paste one of these URLs:
        http://nimbus1:8888/?token=b2a37a0b893cdbb5949aaad1d27e8fcff80841b0a9ae5b1c
     or http://127.0.0.1:8888/?token=b2a37a0b893cdbb5949aaad1d27e8fcff80841b0a9ae5b1c
```
{: .output}

Now, locate the last line in the output, and copy in the clipboard the alphanumeric string that follows after `token=`.  
Open your web browser, and type as URL `<Singularity machine IP Address>:8888`. Paste the clipboard content in the *Token* field and click on *Login*.  

We're in the Jupyter Notebook!  
From the list of files & directories, click on `LogisticRegression.ipynb`, and we'll go through the notebook.

When we're done, press twice `Ctrl-C` to kill the Notebook.


### Running Python code without Jupyter

Because our Jupyter container has a functioning version of Python installed, we can also use it to run our Python code directly from the command line:

```
$ singularity exec plotly.sif python ./logreg.py
```

As a result you should now have a plot saved as a `.png` file.


### Bonus: another Machine Learning Example

Let's try another ML example.  We'll build another Docker image that uses Jupyter.  For this  this example, `cd` to the `demos/07_python/image-classification`.

In here we have some image data (`dataset/`), a Jupyter notebook, and a Dockerfile.

Just like before, we'll have you build your own Jupyter container.

> Build a custom Jupyter/Python container
> For this example you should try building your own Docker container.  Some hints:
>
> * Decide if you should use a base image (the [Jupyter Docker Stacks](https://jupyter-docker
stacks.readthedocs.io/en/latest/using/selecting.html) might be a good place to start)
> * You'll need to install 2 modules: **mahotas** and **opencv-python**
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
