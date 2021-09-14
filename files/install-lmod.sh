#!/bin/bash
# following Lmod documentation

LMOD_VERSION="8.5"

# install pre-requisites
sudo apt-get update
sudo apt-get install -y \
  lua5.3 \
  lua-bit32:amd64 \
  lua-posix:amd64 \
  lua-posix-dev \
  liblua5.3-0:amd64 \
  liblua5.3-dev:amd64 \
  tcl \
  tcl-dev \
  tcl8.6 \
  tcl8.6-dev:amd64 \
  libtcl8.6:amd64

# patch for default Lua tools
sudo update-alternatives --install /usr/bin/lua \
  lua-interpreter /usr/bin/lua5.3 130 \
  --slave /usr/share/man/man1/lua.1.gz lua-manual \
  /usr/share/man/man1/lua5.3.1.gz
sudo update-alternatives --install /usr/bin/luac \
  lua-compiler /usr/bin/luac5.3 130 \
  --slave /usr/share/man/man1/luac.1.gz lua-compiler-manual \
  /usr/share/man/man1/luac5.3.1.gz
sudo ln -s /usr/lib/x86_64-linux-gnu/liblua5.3-posix.so \
  /usr/lib/x86_64-linux-gnu/lua/5.3/posix.so

# install Lmod
wget https://sourceforge.net/projects/lmod/files/Lmod-${LMOD_VERSION}.tar.bz2
tar xf Lmod-${LMOD_VERSION}.tar.bz2
cd Lmod-${LMOD_VERSION}/
./configure --prefix=/opt/apps
sudo make install
cd ..

# configure Lmod
sudo ln -s /opt/apps/lmod/lmod/init/profile /etc/profile.d/z00_lmod.sh
. /etc/profile.d/z00_lmod.sh

