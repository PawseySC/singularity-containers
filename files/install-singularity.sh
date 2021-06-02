#!/bin/bash

# define custom variables
INSTALL_DIR="/usr/local"
GO_VERSION="1.14.12"
SING_VERSION="3.7.4" # or another tag or branch if you like

USERID="$USER" # DO NOT change this
cd $HOME

# install pre-requisites
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup

# install go
OS="linux"
ARCH="amd64"
if [ ! -e go$GO_VERSION.$OS-$ARCH.tar.gz ] ; then
  wget https://dl.google.com/go/go$GO_VERSION.$OS-$ARCH.tar.gz
fi
sudo tar -C $INSTALL_DIR -xzvf go$GO_VERSION.$OS-$ARCH.tar.gz
export PATH="/usr/local/go/bin:${PATH}"
echo "export PATH=/usr/local/go/bin:\${PATH}" >> $(eval echo ~${USERID})/.bashrc
# rm go$GO_VERSION.$OS-$ARCH.tar.gz

# install singularity
if [ ! -e singularity-$SING_VERSION.tar.gz ] ; then
  wget https://github.com/hpcng/singularity/releases/download/v$SING_VERSION/singularity-$SING_VERSION.tar.gz
fi
tar -xzf singularity-$SING_VERSION.tar.gz
cd singularity
./mconfig --prefix=$INSTALL_DIR && \
    make -C ./builddir && \
    sudo make -C ./builddir install
cd ..
rm -r singularity # singularity-$SING_VERSION.tar.gz


# configure singularity
sudo sed -i 's/^ *mount *home *=.*/mount home = no/g' $INSTALL_DIR/etc/singularity/singularity.conf

echo ". ${INSTALL_DIR}/etc/bash_completion.d/singularity" >> $(eval echo ~${USERID})/.bashrc
