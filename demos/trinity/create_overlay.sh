#!/bin/bash

export COUNT="200"
export BS="1M"
export FILE="my_overlay"
singularity exec docker://ubuntu:18.04 bash -c " \
    mkdir -p overlay_tmp/upper && \
    dd if=/dev/zero of=$FILE count=$COUNT bs=$BS && \
    mkfs.ext3 -d overlay_tmp $FILE && \
    rm -rf overlay_tmp \
    "
