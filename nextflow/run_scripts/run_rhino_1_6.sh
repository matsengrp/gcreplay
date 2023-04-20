#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH
nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest_1_6.csv" \
        --results "results/$(date -I)-PR1.6-rhino" \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -resume
