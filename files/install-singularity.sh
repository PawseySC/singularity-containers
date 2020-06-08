#!/bin/bash

# define custom variables
GO_DIR="/opt"
GO_VERSION="1.13.10"

SING_DIR="/usr/local"
SING_VERSION="v3.5.2" # or another tag or branch if you like

USERID="$USER" # do not change this


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


# install singularity
go get -d github.com/sylabs/singularity
cd $GOPATH/src/github.com/sylabs/singularity && \
    git fetch && \
    git checkout $SING_VERSION
./mconfig --prefix=$SING_DIR && \
    make -C ./builddir && \
    sudo make -C ./builddir install


# configure singularity
sudo sed -i 's/^ *mount *home *=.*/mount home = no/g' $SING_DIR/etc/singularity/singularity.conf

echo ". ${SING_DIR}/etc/bash_completion.d/singularity" >> $(eval echo ~${USERID})/.bashrc
