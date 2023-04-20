#!/bin/bash

set -e

nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest_1_6.csv" \
        --results "results/$(date -I)-test" \
        -work-dir "work/$(date -I)-test" \
        -profile docker \
        -resume

#--dms_vscores "${PWD}/data/_ignore/dms/DMS_final_variant_scores.csv" \
#--dms_sites "${PWD}/data/_ignore/dms/DMS_sites.csv" \
