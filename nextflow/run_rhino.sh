#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

# --manifest "data/_ignore/miseq/manifest.csv" \
nextflow run main.nf \
        --plate_barcodes "data/barcodes/plateBC.txt" \
        --well_barcodes "data/barcodes/96FBC.txt" \
        --results "results/$(date -I)" \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-dag "results/$(date -I)/dag.svg" \
        -resume
