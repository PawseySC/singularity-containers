---
title: "Setup Singularity on your machine"
teaching: 15
exercises: 0
questions:
objectives:
- Get an overview of how to install Singularity on various operating systems
keypoints:
- bb
---


### Note on the installation

We're here briefly outlining how to install Singularity on Linux, macOS and Windows.  

In all cases, you will need *admin* privileges in the machine to finalise the installation.


### Linux installation: template script

Let's go through stages to comment on the script that ships with this tutorial: [install-singularity.sh]({{ page.root }}/files/install-singularity.sh).  The script is suitable for Ubuntu machines, and will need modifications for other distributions, especially to install some pre-requisite tools.

The script has been mostly inspired by the [Sylabs docs for installing Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps).

The first part of the script has some variable definitions, some of which might need customisation, *i.e.* version and installation path for Singularity and Go (a build pre-requisite):

```
#!/bin/bash

# define custom variables
GO_DIR="/opt"
GO_VERSION="1.13.10"

SING_DIR="/usr/local"
SING_VERSION="v3.5.2" # or another tag or branch if you like

USERID="$USER" # do not change this
```
{: .source}

Then, we're using the Ubuntu package manager, `apt`, to install some pre-requisite tools.  If running on a different Linux flavour, you might need to modify this part.

```
# install pre-requisites
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    git \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config
```
{: .source}

Next, we're installing *Go*, the build tool that will be used to actually compile the Singularity source code into an executable:

```
# install go
GOPATH="${GO_DIR}/go"
sudo mkdir -p $GOPATH
sudo chown ${USERID}:${USERID} $GOPATH
cd $GO_DIR
OS="linux"
ARCH="amd64"
sudo wget https://dl.google.com/go/go$GO_VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$GO_VERSION.$OS-$ARCH.tar.gz
PATH="/usr/local/go/bin:${PATH}:${GOPATH}/bin"
echo "export GOPATH=$(pwd)/go" >> $(eval echo ~${USERID})/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> $(eval echo ~${USERID})/.bashrc
go get -u github.com/golang/dep/cmd/dep
```
{: .source}

It's now time to download and build Singularity:

```
# install singularity
go get -d github.com/sylabs/singularity
cd $GOPATH/src/github.com/sylabs/singularity && \
    git fetch && \
    git checkout $SING_VERSION
./mconfig --prefix=$SING_DIR && \
    make -C ./builddir && \
    sudo make -C ./builddir install
```
{: .source}

We're almost there!  The last step is about customising the installation:

```
# configure singularity
sudo sed -i 's/^ *mount *home *=.*/mount home = no/g' $SING_DIR/etc/singularity/singularity.conf

echo ". ${SING_DIR}/etc/bash_completion.d/singularity" >> $(eval echo ~${USERID})/.bashrc
```
{: .source}

In this case, we're only doing two things:
* disabling bind mounting of `$HOME` in the configuration file, to enforce best practices;
* enabling bash shell completion, to make our life easier when using Singularity.

The Singularity installation is highly customisable through the configuration file `$SING_DIR/etc/singularity/singularity.conf`.  Refer to the [Sylabs documentation](https://sylabs.io/guides/3.5/user-guide/index.html) for more information.


### macOS and Windows installation, option 1: using Vagrant

