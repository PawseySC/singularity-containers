#!/bin/bash

export PASSWORD=rstudiopassword
echo $USER && echo $PASSWORD

singularity instance start \
    -c \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    tidyverse_long.sif \
    myserver