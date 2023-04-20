#!/bin/bash

set -e

source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH
module load git-lfs/2.11.0

nextflow run smartseq.nf \
        --hk_df data/smartseq/2022-07-10-smartseq-plus-miseq-hk.csv \
        --results "results/2022-07-11-smartseq-miseq" \
        -profile fred_hutch_rhino \
        -work-dir /fh/scratch/delete30/matsen_e/jgallowa/work/ \
        -with-trace \
        -resume

git add results/2022-07-11-smartseq-miseq/merged-results/*.csv
git add results/2022-07-11-smartseq-miseq/gctrees/PR1.*GC/gctree.p
git add results/2022-07-11-smartseq-miseq/gctrees/PR1.*GC/L*.svg
git add results/2022-07-11-smartseq-miseq/gctrees/PR1.*GC/delta*.svg
git add results/2022-07-11-smartseq-miseq/gctrees/PR1.*GC/*.nk

