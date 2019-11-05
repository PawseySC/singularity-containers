#!/bin/bash

MPICH_VERSION="3.1.4"
MPICH_CONFIGURE_OPTIONS="--enable-fast=all,O3 --prefix=$(pwd)/mpich-${MPICH_VERSION}/apps"
MPICH_MAKE_OPTIONS="-j4"
unset F90

wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz
tar xvzf mpich-${MPICH_VERSION}.tar.gz

cd mpich-${MPICH_VERSION}

./configure ${MPICH_CONFIGURE_OPTIONS}
make ${MPICH_MAKE_OPTIONS}
make install

