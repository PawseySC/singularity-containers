#!/bin/bash

enroot_ver="3.0.2"
enroot_flavour="-hardened"

sudo apt-get update
sudo apt-get install -y curl gawk jq squashfs-tools parallel

curl -fSsL -O https://github.com/NVIDIA/enroot/releases/download/v${enroot_ver}/enroot${enroot_flavour}_${enroot_ver}-1_amd64.deb
curl -fSsL -O https://github.com/NVIDIA/enroot/releases/download/v${enroot_ver}/enroot${enroot_flavour}+caps_${enroot_ver}-1_amd64.deb # optional
chmod o+r enroot${enroot_flavour}*.deb
sudo apt-get install -y ./enroot${enroot_flavour}*.deb
