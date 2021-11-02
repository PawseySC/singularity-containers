#!/bin/bash

# essential for reproducibility of installation - leave blank "" to use latest commit
shpc_checkout="0.0.34"

# some pre-requisites are assumed available: singularity and lmod - see their install scripts
#
# install other pre-requisites
sudo apt-get update
sudo apt-get install -y \
  python3 \
  python3-pip \
  python3-setuptools
pip3 install --upgrade pip

# internal variables - do not edit
USERID="$USER"
python_ver="$( python3 -V | cut -d ' ' -f 2 | cut -d . -f 1,2 )"
install_dir="$(pwd)/singularity-hpc"
lib_dir="$install_dir/lib/python${python_ver}/site-packages"
bin_dir="$install_dir/bin"
mod_dir="$install_dir/modules"
#
export PYTHONPATH="$lib_dir:$PYTHONPATH" 
export PATH="$bin_dir:$PATH"

# install shpc
git clone https://github.com/singularityhub/singularity-hpc
cd singularity-hpc
git checkout $shpc_checkout
pip3 install -e . --prefix="$(pwd)"
cd ..

# configure shpc
mkdir -p $mod_dir
module use $mod_dir
shpc config set container_base:\$root_dir/containers
echo "export PYTHONPATH=\"$lib_dir:\$PYTHONPATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "export PATH=\"$bin_dir:\$PATH\"" >> $(eval echo ~${USERID})/.bashrc
echo "module use $mod_dir" >> $(eval echo ~${USERID})/.bashrc

