#!/bin/bash
# use docker://trinityrnaseq/trinityrnaseq:2.8.6

Trinity \
    --seqType fq --left trinity_test_data/reads.left.fq.gz  \
    --right trinity_test_data/reads.right.fq.gz \
    --max_memory 1G --CPU 1 --output <OUTPUT-DIRECTORY>
