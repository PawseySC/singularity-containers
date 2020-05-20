#!/bin/bash

USERID="ubuntu"

CHARLIE_VER="0.15"

# install Spack dependencies
sudo apt update
sudo apt install -y \
  python3 \
  build-essential \
  make \
  git \
  curl

# create install dir
export SPACK_ROOT="/opt/spack"
sudo mkdir -p $SPACK_ROOT
sudo chown ${USERID}:${USERID} $SPACK_ROOT

# clone Spack
cd $SPACK_ROOT
git clone https://github.com/spack/spack.git .

# configure shell environment for Spack
echo ". ${SPACK_ROOT}/share/spack/setup-env.sh" >> ~/.bashrc
. ${SPACK_ROOT}/share/spack/setup-env.sh


# now let's install Charliecloud through Spack
# this will take a while ..
spack install charliecloud@${CHARLIE_VER}

