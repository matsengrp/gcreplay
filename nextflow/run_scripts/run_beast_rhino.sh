#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity

export PATH=$SINGULARITYROOT/bin/:$PATH

results="results/$(date -I)-beast"
nextflow run beast.nf \
    --seqs "results/gcreplay-observed/*-*-*-*-[1-6]-*-*-*.fasta" \
    --chain_length 50000000 \
    --log_every 10000 \
    -profile fred_hutch_rhino \
    --burn_frac 0.0 \
    --results $results \
    -work-dir $results/work \
    -with-trace $results/trace.txt \
    -with-report $results/report.html \
    --save_pkl_trees false \
    -resume
    # --chain_length 100000 \
    # --log_every 1000 \

