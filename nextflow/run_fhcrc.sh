#!/bin/bash

set -euo pipefail
source /app/lmod/lmod/init/profile

# Nextflow Version
NXF_VER=21.04.3

module load Nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH
results_dir="results/"
workdir="/fh/scratch/delete30/matsen_e/jgallowa/work"
manifest="manifest.csv"
prefix="data/input"

NXF_VER=$NXF_VER \
nextflow run main.nf \
        --manifest $manifest \
        --reads_prefix $prefix \
        --results $results_dir \
        --dms_vscores "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/improved-Kd-fitting/tite-seq-modeling/output/final_variant_scores.csv" \
        -work-dir $workdir \
        -profile fred_hutch_rhino \
        -resume
