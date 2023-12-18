#!/bin/bash

set -euo pipefail
# source /app/lmod/lmod/init/profile

# Nextflow Version
NXF_VER=21.04.3

module load Nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH
results_dir="results/$(date -I)-PR2.1-rhino"
workdir="/fh/scratch/delete30/matsen_e/jgallowa/work"
manifest="data/_ignore/Early_TP_Single_GC_PR_2/manifest_PR2.csv"

NXF_VER=$NXF_VER \
nextflow run main.nf \
        --manifest $manifest \
        --results $results_dir \
        -work-dir $workdir \
        -profile fred_hutch_rhino \
        -resume