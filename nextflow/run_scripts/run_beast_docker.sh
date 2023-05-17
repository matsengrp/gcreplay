#!/bin/bash

set -e

nextflow run beast.nf \
    -with-report \
    -profile docker \
    --results results/$(date -I)-beast \
    --chain_length 100000 \
    --log_every 1000 \
    --burn_frac 0.0 \
    -resume
