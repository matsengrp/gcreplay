#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH
for run in 3 4;
do
nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest_1_6.csv" \
        --results "results/$(date -I)P1.6-run${run}-seed23" \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-dag "results/$(date -I)/dag.svg"
done
