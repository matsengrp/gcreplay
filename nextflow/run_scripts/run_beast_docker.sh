#!/bin/bash

set -e

nextflow run beast.nf \
    -with-report \
    -profile docker \
    --results results/$(date -I)-beast-alloc-test \
    --chain_length 100000 \
    --log_every 1000 \
    --burn_frac 0.0 \
    --save_pkl_trees false
    # -resume
