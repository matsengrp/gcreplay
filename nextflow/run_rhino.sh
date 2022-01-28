#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow run main.nf \
        --manifest data/_ignore/miseq/manifest.csv \
        --results "$(date -I)" \
        -profile fred_hutch_rhino \
        -with-dag "$(date -I)/dag.svg" \
        -resume
