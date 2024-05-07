#!/bin/bash

set -euo pipefail
# source /app/lmod/lmod/init/profile

# Nextflow Version
#NXF_VER=21.04.3

# Load the Nextflow module (if running on rhino/gizmo)
ml nextflow
# Load the Apptainer module (if running on rhino/gizmo with Apptainer)
ml Apptainer
export PATH=$APPTAINERROOT/bin/:$PATH
results_dir="results"
workdir="/fh/scratch/delete30/matsen_e/jgallowa/work"
manifest="manifest.csv"
prefix="data/input"

# export NXF_VER=21.04.3
# echo $NXF_VER
nextflow run main.nf \
        --manifest $manifest \
        --reads_prefix $prefix \
        --results $results_dir \
        -work-dir $workdir \
        -profile fred_hutch_rhino \
        -resume
