#!/bin/bash

set -e

nextflow run simulations.nf \
        --results "results/$(date -I)-sim-docker" \
        --test true \
        -work-dir work-sim \
        -profile docker \
        -with-dag "results/$(date -I)/sim-dag.svg" \
        -resume

#--dms_vscores "${PWD}/data/_ignore/dms/DMS_final_variant_scores.csv" \
#--dms_sites "${PWD}/data/_ignore/dms/DMS_sites.csv" \

