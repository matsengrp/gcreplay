#!/bin/bash

set -e

nextflow run smartseq.nf \
        --results "results/$(date -I)-smartseq-docker" \
        --dms_vscores "${PWD}/data/_ignore/dms/DMS_final_variant_scores.csv" \
        --dms_sites "${PWD}/data/_ignore/dms/DMS_sites.csv" \
        -work-dir work-ss \
        -profile docker \
        -with-dag "results/$(date -I)/ss-dag.svg" \
        -resume

