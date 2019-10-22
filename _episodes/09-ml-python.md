---
title: "Making Python not awful with containers"
teaching: 0
exercises: 20
questions:
objectives:
- Embed a Python app in a container and run it
keypoints:
- Containers are great way to manage Python workflows
---

### Why can Python be Awful? ###

Python is a great language for doing all kinds of work.

![Python can get messy]({{ page.root }}/fig/python-complexity-cartoon.png)


### Build and run a dockerised Python app ###

This is a quick Python data science example.  We'll build a Python container, add a simple Python app, and run it.

To begin, let's clone another repo:

```
$ git clone https://github.com/skjerven/python-demo.git
$ cd python-demo
```
{: .bash}

There a few files here:

* `Dockerfile`: outlines how we'll build our container

```
# Use an official Python runtime as a parent image
FROM python:3.6-slim

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Define environment variable
ENV NAME World

# Run app.py when the container launches
CMD ["python", "my_app.py"]
```
{: .source}

* `requirements.txt`: what python packages we want to install with pip

```
numpy
scipy
scikit-learn
```
{: .source}

* `my_app.py`: a simple python app that builds a decision tree using Scikit-Learn

```
from sklearn import tree
#DataSet
#[size,weight,texture]
X = [[181, 80, 44], [177, 70, 43], [160, 60, 38], [154, 54, 37],[166, 65, 40],
     [190, 90, 47], [175, 64, 39],
     [177, 70, 40], [159, 55, 37], [171, 75, 42], [181, 85, 43]]

Y = ['apple', 'apple', 'orange', 'orange', 'apple', 'apple', 'orange', 'orange',
     'orange', 'apple', 'apple']

#classifier - DecisionTreeClassifier
clf_tree = tree.DecisionTreeClassifier();
clf_tree = clf_tree.fit(X,Y);

#test_data
test_data = [[190,70,42],[172,64,39],[182,80,42]];

#prediction
prediction_tree = clf_tree.predict(test_data);

print("Prediction of DecisionTreeClassifier:",prediction_tree);
```
{: .python}

There are some aspects of this Dockerfile that are worth mentioning:

* we are using the base image `python:3.6-slim`, which comes from the official Python repository on Docker Hub; images in this repo whose tag has `slim` are characterised by a very small size, in this case about 100 MB as opposed to about 400 MB for the `miniconda3` image;
* the Docker instruction `ADD` permits to copy the content of the build context from the host machine inside the container image;
* the `CMD` instruction sets the default behaviour for this container: running the embedded Python app;
* finally, note from the requirement libraries that this is quite a standard scientific Python stack.

Now to build our container we can simply run:

```
$ docker build -t python-demo .
```
{: .bash}

After that, we can run it with:

```
$ docker run python-demo
```
{: .bash}


### Run a Python app on HPC with Shifter ###

#### a) re-use the image we have just built ####

On our Docker machine let us push the container image we created above to Docker Hub; you'll need a free account.

First we must give the image a name that complies with the Hub's nomenclature (see previous episode on build):

```
$ docker tag python-demo <your-dockerhub-account>/python-demo
```
{: .bash}

Now let's push the image:

```
$ docker push <your-dockerhub-account>/python-demo
```
{: .bash}

```
The push refers to repository [docker.io/marcodelapierre/python-demo]
862d6710cd15: Pushed 
302ce4960403: Pushed 
[..]
latest: digest: sha256:4db5f0f69cc888d47f4c4b4cac33fad6b004a8e333b36a699ebd43f5b44a7241 size: 1999
```
{: .output}

We are now moving to the Pawsey HPC system. Let's pull the image, then change directory to either `$MYSCRATCH` or `$MYGROUP`:

```
$ module load shifter
$ sg $PAWSEY_PROJECT -c 'shifter pull <your-dockerhub-account>/python-demo'

$ cd $MYSCRATCH
```
{: .bash}

Let us write a SLURM script to execute our Python app using this container, we'll use our favourite text editor to save it as `python.sh` (remember to specify your Pawsey project ID in the script!): 

```
#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --export=NONE
#SBATCH --job-name=python

module load shifter

# run Python app
srun --export=all shifter run <your-dockerhub-account>/python-demo
```
{: .bash}

Let's submit it to the SLURM scheduler:

```
$ sbatch --reservation <your-pawsey-reservation> python.sh
```
{: .bash}

#### b) use a publicly available image for scientific Python ####

In this case we are going to use `jupyter/scipy-notebook`:

```
$ module load shifter
$ sg $PAWSEY_PROJECT -c 'shifter pull jupyter/scipy-notebook'

$ cd $MYSCRATCH
```
{: .bash}

Let us write a second SLURM script, we'll call it `python2.sh`. Contrary to Docker example above, our Python app is not embedded in the image, so we'll need to explicitly download it from the Git repo, and then run it through the Python interpreter in the container:

```
#!/bin/bash -l

#SBATCH --account=<your-pawsey-project>
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --export=NONE
#SBATCH --job-name=python2

module load shifter

# clone Git repo with the app
git clone https://github.com/skjerven/python-demo.git
cd python-demo

# run Python app
srun --export=all shifter run jupyter/scipy-notebook python my_app.py
```
{: .bash}

Finally we are submitting the script with SLURM:

```
$ sbatch --reservation <your-pawsey-reservation> python2.sh
```
{: .bash}

