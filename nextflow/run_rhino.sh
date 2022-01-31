#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

#--manifest data/_ignore/miseq/manifest.csv \
nextflow run main.nf \
        --reads "data/_ignore/miseq/PR1.6/*_R{1,2}*.fastq" \
        --results "$(date -I)" \
        -profile fred_hutch_rhino \
        -with-dag "$(date -I)/dag.svg" \
        -resume
