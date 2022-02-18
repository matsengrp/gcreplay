#!/bin/bash
set -eu

# ------- 
# Post annotation wrangle. This comed directly after the
# partis annotation

# 

IGH_AIRR=2022-02-09/partis_annotation/PR-1-6/engrd/single-chain/partition-igh.tsv
IGK_AIRR=2022-02-09/partis_annotation/PR-1-6/engrd/single-chain/partition-igk.tsv

# TODO pass this down in the pipeline
FASTA=2022-02-09/ranked_bcr_sequences_per_well/PR-1-6.fasta

# TODO change the key that's getting passed down with 
# with the pipeline currently to the key_file?
KEY_FILE=PR1_6-key.csv
OUT=gc-df.csv

# wrangle annotation -> gc merged dataframe
python gcreplay-tools.py wrangle-annotation \
    --igh-airr $IGH_AIRR \
    --igk-airr $IGK_AIRR \
    --input-fasta $FASTA \
    --key-file $KEY_FILE \
    -o $OUT

echo Done with wrangle, writing df to $OUT
