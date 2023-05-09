#!/bin/bash

set -e

nextflow run beast.nf \
    -with-report \
    -profile docker \
    --results results/$(date -I)-beast \
    --chain_length 10000 \
    --log_every 1000 \
    -resume
