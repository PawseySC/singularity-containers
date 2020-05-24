#!/bin/bash

export PASSWORD=rstudiopassword
echo $USER && echo $PASSWORD

singularity exec \
    -c \
    -B $(pwd):/home/rstudio \
    -B $(pwd):$HOME \
    tidyverse_3.6.1.sif \
    rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=0 --auth-pam-helper-path=pam-helper