#!/bin/bash

export COUNT="200"
export BS="1M"
singularity exec docker://ubuntu:18.04 bash -c " \
    mkdir -p overlay_tmp/upper && \
    dd if=/dev/zero of=my_overlay count=$COUNT bs=$BS && \
    mkfs.ext3 -d overlay_tmp my_overlay && \
    rm -rf overlay_tmp \
    "
