#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --results "results/$(date -I)" \
        --bcr_count_thresh 5 \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-dag "results/$(date -I)/dag.svg" \
        -with-trace \
        -resume
