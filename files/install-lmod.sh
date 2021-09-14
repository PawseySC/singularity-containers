#!/bin/bash

LMOD_VER="8.3"
SPACK_ROOT="/opt/spack"

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

# create install dir
sudo mkdir -p $SPACK_ROOT
sudo chown ${USERID}:${USERID} $SPACK_ROOT

# clone Spack
cd $SPACK_ROOT
git clone https://github.com/spack/spack.git .

# configure shell environment for Spack
echo ". ${SPACK_ROOT}/share/spack/setup-env.sh" >> $(eval echo ~${USERID})/.bashrc
. ${SPACK_ROOT}/share/spack/setup-env.sh


# now let's install Lmod through Spack
# this will take a while ..
spack install lmod@${LMOD_VER}
LMOD_DIR="$(spack location -i lmod)"
echo ". ${LMOD_DIR}/lmod/lmod/init/profile" >> $(eval echo ~${USERID})/.bashrc
. ${LMOD_DIR}/lmod/lmod/init/profile
