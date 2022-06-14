#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
module load git-lfs/2.11.0

export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow run main.nf \
        --manifest "data/_ignore/miseq/manifest.csv" \
        --results "results/$(date -I)" \
        --bcr_count_thresh 5 \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -profile fred_hutch_rhino \
        -with-dag "results/$(date -I)/dag.svg" \
        -with-trace \
        -resume

rm results/latest
ln -s results/$(date -I) results/latest

git add results/$(date -I)/merged-results/*.csv
git add results/$(date -I)/gctrees/PR1.*GC/gctree.p
git add results/$(date -I)/gctrees/PR1.*GC/L*.svg
git add results/$(date -I)/gctrees/PR1.*GC/delta*.svg
git add results/$(date -I)/gctrees/PR1.*GC/*.nk


