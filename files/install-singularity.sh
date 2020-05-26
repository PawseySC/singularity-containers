#!/bin/bash

export USERID=$USER

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
export GO_DIR=/opt
export GOPATH=${GO_DIR}/go
sudo mkdir -p $GOPATH
sudo chown ${USERID}:${USERID} $GOPATH
cd $GO_DIR
export VERSION=1.13.10 OS=linux ARCH=amd64
sudo wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
export GOPATH=$(pwd)/go
export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin
echo "export GOPATH=$(pwd)/go" >> $(eval echo ~${USERID})/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> $(eval echo ~${USERID})/.bashrc
go get -u github.com/golang/dep/cmd/dep

# install singularity
export SING_DIR=/usr/local
export SING_VERSION=v3.5.2 # or another tag or branch if you like && \
go get -d github.com/sylabs/singularity
cd $GOPATH/src/github.com/sylabs/singularity && \
    git fetch && \
    git checkout $SING_VERSION
./mconfig --prefix=$SING_DIR && \
    make -C ./builddir && \
    sudo make -C ./builddir install
sudo sed -i 's/^ *mount *home *=.*/mount home = no/g' $SING_DIR/etc/singularity/singularity.conf

echo ". ${SING_DIR}/etc/bash_completion.d/singularity" >> $(eval echo ~${USERID})/.bashrc

