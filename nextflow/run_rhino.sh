#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
module load git-lfs/2.11.0

export PATH=$SINGULARITYROOT/bin/:$PATH
results="results/$(date -I)"
nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --results $results \
        --bcr_count_thresh 5 \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-trace $results/trace.txt \
        -with-report $results/report.html \
        -resume

(cd results && rm latest && ln -s $(date -I) latest && git add latest)

git add results/$(date -I)/merged-results/*.csv
git add results/$(date -I)/gctrees/PR1.*GC/gctree.p
git add results/$(date -I)/gctrees/PR1.*GC/L*.svg
git add results/$(date -I)/gctrees/PR1.*GC/delta*.svg
git add results/$(date -I)/gctrees/PR1.*GC/*.nk
