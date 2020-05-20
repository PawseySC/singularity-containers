#!/bin/bash

sudo apt update
sudo apt install -y environment-modules

cd /usr/bin

# assuming version 8.6 here
sudo ln -s tclsh8.6 tclsh

cd -