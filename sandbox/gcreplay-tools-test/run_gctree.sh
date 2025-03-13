#! /usr/bin/env bash
set -eu


# -------
# There is a split csv and the gc is sent combined with full gc df
# it then queries the df, makes a fasta, runs gctree, and finally does a second wrangle
# THIS IS IF GC != 0
# The other case will 

GCDF=gc-df.csv                  # Full PR cell HK database (csv)
GC_NUM_HC=87                    # GC you would like to extract
QUERY="GC_num_HC==${GC_NUM_HC}" # The query on the dataset before fasta extraction

# query the gc dataframe
python gcreplay-tools.py query-df -df $GCDF -q $QUERY -o gc_${GC_NUM_HC}_hk.csv
echo done: query

# convert some columns into headers and sequences for a basic fasta format conversion
python gcreplay-tools.py gc-df-to-fasta --gc-df gc_${GC_NUM_HC}_hk.csv -o gc_${GC_NUM_HC}_hk.fasta
echo done: created fasta

## First we deduplicate the sequences and convert to phylip alignment format, 
## and also determine sequence abundances. The deduplicate command writes 
## the phylip alignment of unique sequences to stdout (which we redirect to a file here). 
deduplicate gc_${GC_NUM_HC}_hk.fasta \
    --root naive \
    --idmapfile idmap-gc${GC_NUM_HC}-HK.txt \
    --abundance_file abundances.csv > deduplicated.phylip
echo done: deduplicate

# TODO describe
mkconfig deduplicated.phylip dnapars > dnapars.cfg
echo done: mkconfig

# TODO describe
rm outtree outfile
dnapars < dnapars.cfg > dnapars.log
echo done: dnapars

# TODO describe
gctree infer outfile abundances.csv \
    --root naive \
    --frame 1 \
    --frame2 1 \
    --positionmapfile ../../data/position_maps/210929PR_H_index-imgtTs.txt \
    --positionmapfile2  ../../data/position_maps/210930PR_K_index-imgtTs.txt \
    --chain_split 336 \
    --idlabel \
    --verbose \
    | tee gctree.inference.log
echo done: gctree

# TODO describe
## Post-GCtree processing
#- import idmap, add to GCdf
#- import putative nodes, fill idmap info, add to GCdf
#- calculate total mut#
#- calculate DMS scores
