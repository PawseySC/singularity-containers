#!/bin/bash

export SINGULARITYENV_USER=$USER
export SINGULARITYENV_PASSWORD=rstudiopassword
echo $SINGULARITYENV_USER && echo $SINGULARITYENV_PASSWORD

singularity instance start \
    -C \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    tidyverse_long.sif \
    myserver