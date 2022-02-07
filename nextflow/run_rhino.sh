#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --plate_barcodes "data/barcodes/plateBC.txt" \
        --well_barcodes "data/barcodes/96FBC.txt" \
        --results "results/$(date -I)" \
        -profile fred_hutch_rhino \
        -with-dag "$(date -I)/dag.svg" \
        -resume
