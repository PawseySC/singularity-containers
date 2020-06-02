#!/bin/bash

export USERID=$USER

sudo apt update
sudo apt install -y environment-modules

cd /usr/bin
# assuming version 8.6 here
sudo ln -s tclsh8.6 tclsh

MODULE_VERSION=$( module -V | cut -d ' ' -f 3 )
echo "export MODULE_VERSION=\"$MODULE_VERSION\"" >> $(eval echo ~${USERID})/.bashrc

