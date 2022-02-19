#!/bin/bash
set -eu

# ------- 
# Post annotation wrangle. This comed directly after the
# partis annotation

# 

for PR in 1 2 3 4 5 6 7 8
do
    # PR=6
    RESDIR=../../results/2022-02-09/
    IGH_AIRR=$RESDIR/partis_annotation/PR-1-$PR/engrd/single-chain/partition-igh.tsv
    IGK_AIRR=$RESDIR/partis_annotation/PR-1-$PR/engrd/single-chain/partition-igk.tsv

    # TODO pass this down in the pipeline
    FASTA=$RESDIR/ranked_bcr_sequences_per_well/PR-1-$PR.fasta

    # TODO change the key that's getting passed down with 
    # with the pipeline currently to the key_file?
    KEY_FILE=../../data/key_files/PR1_${PR}-key.csv
    OUT=gc_dfs/gc-df-1-$PR.csv

    # wrangle annotation -> gc merged dataframe
    python gcreplay-tools.py wrangle-annotation \
        --igh-airr $IGH_AIRR \
        --igk-airr $IGK_AIRR \
        --input-fasta $FASTA \
        --key-file $KEY_FILE \
        -o $OUT

    echo Done with wrangle, writing df to $OUT
done
