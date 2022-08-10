#!/bin/bash

set -e

nextflow run simulations.nf \
        --simulations "data/simulations/*.fasta" \
        --results "results/$(date -I)-sim-docker" \
        -work-dir work-sim \
        -profile docker \
        -with-dag "results/$(date -I)/sim-dag.svg"
