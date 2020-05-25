#!/bin/bash

export SINGULARITYENV_USER=$USER
export SINGULARITYENV_PASSWORD=rstudiopassword
echo $SINGULARITYENV_USER && echo $SINGULARITYENV_PASSWORD

singularity exec \
    -C \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    docker://rocker/tidyverse:3.6.1 \
    rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=0 --auth-pam-helper-path=pam-helper
