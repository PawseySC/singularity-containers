#!/bin/bash

NEXTFLOW_VER="19.10.0"

sudo apt-get update
sudo apt-get install -y default-jre

NEXTFLOW_DIR="/usr/local/bin"
cd $NEXTFLOW_DIR
sudo wget https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VER}/nextflow
sudo chmod a+rx nextflow
cd ~
${NEXTFLOW_DIR}/nextflow info

