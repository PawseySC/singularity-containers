#!/bin/bash

MPICH_VERSION="3.1.4"

USERID="$USER"

MPICH_ROOT="/opt/mpich"
MPICH_DIR="$MPICH_ROOT/mpich-$MPICH_VERSION/apps"
MPICH_CONFIGURE_OPTIONS="--enable-fast=all,O3 --prefix=$MPICH_DIR"
MPICH_MAKE_OPTIONS="-j4"
unset F90

sudo apt-get update
sudo apt-get install -y \
    build-essential \
    gfortran

sudo mkdir -p $MPICH_ROOT
sudo chown ${USERID}:${USERID} $MPICH_ROOT

cd $MPICH_ROOT

wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz
tar xvzf mpich-${MPICH_VERSION}.tar.gz

cd mpich-${MPICH_VERSION}

./configure ${MPICH_CONFIGURE_OPTIONS}
make ${MPICH_MAKE_OPTIONS}
make install

echo "export PATH=\"$MPICH_DIR/bin:\$PATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "export CPATH=\"$MPICH_DIR/include:\$CPATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "export LD_LIBRARY_PATH=\"$MPICH_DIR/lib:\$LD_LIBRARY_PATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "export LIBRARY_PATH=\"$MPICH_DIR/lib:\$LIBRARY_PATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "export MANPATH=\"$MPICH_DIR/share/man:\$MANPATH\"" >> $(eval echo ~${USERID})/.bashrc

