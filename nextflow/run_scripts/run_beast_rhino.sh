#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity

export PATH=$SINGULARITYROOT/bin/:$PATH

results="results/$(date -I)-beast-1-6"
nextflow run beast.nf \
    --seqs "results/gcreplay-observed/*PR-1-6*.fasta" \
    --chain_length 50000000 \
    --log_every 10000 \
    --burn_frac 0.0 \
    --results $results \
    -work-dir $results/work \
    -profile fred_hutch_rhino \
    -with-trace $results/trace.txt \
    -with-report $results/report.html
    # -resume

