#!/bin/bash

# define variables
USERID="$USER" # do not change this

# remove previous docker installation, if existing
sudo apt-get -y purge docker docker-engine docker.io

# install pre-requisites
sudo apt-get update
sudo apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common \
    wget

# docker repository: get GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
# docker repository: add
sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update

# install docker
sudo apt-get install -y docker-ce

# ensure docker service is started
sudo systemctl start docker
sudo systemctl enable docker

# add user to docker group, to ditch need for explicit sudo
sudo usermod -aG docker $USERID  ## requires logout, login to take effect
