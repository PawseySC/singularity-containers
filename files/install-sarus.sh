#!/bin/bash


# sarus version (versions older than 1.1.0 require a different installation script)
sarus_ver="1.1.0"
sarus_root="/opt/sarus"


# do not change these
USERID="$USER"
sarus_dir="$sarus_root/${sarus_ver}-Release"


# install system prereqs
sudo apt-get update
sudo apt-get install -y squashfs-tools


# install Sarus
sudo mkdir -p $sarus_root
sudo chmod o+rX $sarus_root
cd $sarus_root
sudo wget https://github.com/eth-cscs/sarus/releases/download/${sarus_ver}/sarus-Release.tar.gz
sudo tar xzf sarus-Release.tar.gz
sudo chmod o+rX ${sarus_ver}-Release
cd ${sarus_ver}-Release
sudo ./configure_installation.sh


# add path to ~/.bashrc
echo "export PATH=\"$sarus_dir/bin:\$PATH\"" >> $(eval echo ~${USERID})/.bashrc


# edit configuration (disable home)
sudo sed -i \
  -e 's;\"source\": *\"/home\",;"source": "/data",;' \
  -e 's;\"destination\": *\"/home\",;"destination": "/data",;' \
  $sarus_dir/etc/sarus.json

