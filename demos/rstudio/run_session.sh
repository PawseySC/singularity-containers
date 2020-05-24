#!/bin/bash

export PASSWORD=rstudiopassword
echo $USER && echo $PASSWORD
export HOME_USER=$USER && [ "$(id -u)" == "1000" ] && export HOME_USER=rstudio

singularity exec \
    -c \
    -B $(pwd):/home/$HOME_USER \
    tidyverse_3.6.1.sif \
    rserver --www-port 8787 --www-address 0.0.0.0 --auth-none=0 --auth-pam-helper-path=pam-helper