#!/bin/bash

export PASSWORD=rstudiopassword
echo $USER && echo $PASSWORD
export HOME_USER=$USER && [ "$(id -u)" == "1000" ] && export HOME_USER=rstudio

singularity instance start \
    -c \
    -B $(pwd):/home/$HOME_USER \
    tidyverse_long.sif \
    myserver