#!/bin/bash
# use docker://biocontainers/blast:v2.2.31_cv2

cd ../blast_db
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
cd -

blastp -query P04156.fasta -db $TUTO/demos/blast_db/zebrafish.1.protein.faa -out results.txt
