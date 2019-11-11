#!/bin/bash

mkdir -p $SIFPATH

singularity pull --dir $SIFPATH docker://biocontainers/blast:v2.2.31_cv2
singularity pull --dir $SIFPATH docker://trinityrnaseq/trinityrnaseq:2.8.6
singularity pull --dir $SIFPATH docker://nvcr.io/hpc/gromacs:2018.2
singularity pull --dir $SIFPATH library://marcodelapierre/beta/openfoam:v1812