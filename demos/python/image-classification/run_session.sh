#!/bin/bash

singularity exec \
    -C \
    -B $(pwd):$HOME \
    class.sif \
    jupyter notebook --no-browser --port=8888 --ip 0.0.0.0 --notebook-dir=$HOME
