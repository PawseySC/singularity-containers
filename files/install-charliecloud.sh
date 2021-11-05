#!/bin/bash

CHARLIE_VER="0.25"
SPACK_ROOT="/opt/spack"
SPACK_VER="releases/v0.16"

USERID="$USER"

# install Spack dependencies
sudo apt update
sudo apt install -y \
  python3 \
  build-essential \
  make \
  git \
  curl \
  unzip

# create Spack install dir
sudo mkdir -p $SPACK_ROOT
sudo chown ${USERID}:${USERID} $SPACK_ROOT
# clone Spack
cd $SPACK_ROOT
git clone https://github.com/spack/spack.git .
cd spack
git checkout ${SPACK_VER}
cd ..
# configure Spack shell environment
echo ". ${SPACK_ROOT}/share/spack/setup-env.sh" >> $(eval echo ~${USERID})/.bashrc
. ${SPACK_ROOT}/share/spack/setup-env.sh


# now let's install Charliecloud through Spack
# this will take a while ..
spack install charliecloud@${CHARLIE_VER}

