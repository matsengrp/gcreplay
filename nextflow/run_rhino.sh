#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

#--plate_barcodes "data/barcodes/plateBC.txt" \
#--well_barcodes "data/barcodes/96FBC.txt" \
#--hdag_sub "data/hdag/MK_RS5NF_substitution.csv" \
#--hdag_mut "data/hdag/MK_RS5NF_mutability.csv" \
#--dms_vscores "data/_ignore/dms/DMS_final_variant_scores.csv" \
#--dms_sites "data/_ignore/dms/DMS_sites.csv" \
nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --results "results/$(date -I)" \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-dag "results/$(date -I)/dag.svg" \
        -resume
