#!/bin/bash
set -euo pipefail


# This file is the template for running any particular action
# on the inference, organization with final dataset, node features
# on any one gc. The pipeline will manage the 
# the three parameters below giving you a copy of the full dataset,
# and the gc num that this process is currently running (since it will be in parallel)


echo "Running phipflow"

# should be getting these thins from the pipeline
GCDF=gc-df.csv
SUB=MK_RS5NF_substitution.csv
MUT=MK_RS5NF_mutability.csv
DMS_VSCORES=DMS_final_variant_scores.csv
DMS_SITES=DMS_sites.csv


# will not need to do this in the pipeline 
GC_NUM_HC=87                    # GC you would like to extract
QUERY="GC_num_HC==${GC_NUM_HC}" # The query on the dataset before fasta extraction
python gcreplay-tools.py query-df -df $GCDF -q $QUERY -o gc_${GC_NUM_HC}_hk.csv
echo \(LOG\) done: query


# convert some columns into headers and sequences for a basic fasta format conversion
python gcreplay-tools.py gc-df-to-fasta \
    --gc-df gc_${GC_NUM_HC}_hk.csv \
    -o gc_${GC_NUM_HC}_hk.fasta
echo \(LOG\) done: created fasta


# make isotype map for the specific gc
# TODO how should this be formatted?
python gcreplay-tools.py get-columns \
    -df gc_${GC_NUM_HC}_hk.csv \
    -c ID_HK -c isotype_HC \
    -o gc_${GC_NUM_HC}_hk.isotypemap
echo \(LOG\) done: extracted isotypes
#../gcreplay.isotypemap \

## First we deduplicate the sequences and convert to phylip alignment format, 
## and also determine sequence abundances. The deduplicate command writes 
## the phylip alignment of unique sequences to stdout (which we redirect to a file here). 
deduplicate gc_${GC_NUM_HC}_hk.fasta \
    --root naive \
    --idmapfile idmap-gc${GC_NUM_HC}-HK.txt \
    --abundance_file abundances.csv > deduplicated.phylip
echo \(LOG\) done: deduplicate

# TODO describe
mkconfig deduplicated.phylip dnapars > dnapars.cfg
echo \(LOG\) done: mkconfig

# TODO describe
rm -rf outtree outfile
dnapars < dnapars.cfg > dnapars.log
echo \(LOG\) done: dnapars

# Run inference on sequences, 
mkdir -p gctree_output


xvfb-run -a gctree infer outfile abundances.csv \
    --idmapfile idmap-gc${GC_NUM_HC}-HK.txt \
    --isotype_mapfile ./gc_${GC_NUM_HC}_hk.isotypemap \
    --mutability $MUT \
    --substitution $SUB \
    --root naive \
    --verbose \
    --outbase gctree_output/gctree \
    | tee gctree.inference.log
echo \(LOG\) done: gctree


xvfb-run -a python gcreplay-tools.py featurize-nodes \
    gctree_output/gctree.inference.1.p \
    DMS_final_variant_scores.csv \
    DMS_sites.csv \
    --igk_idx 336 \
    --output_dir gctree_output/
echo \(LOG\) done: Viz

# TODO wrangle in putative nodes and join with gc_${GC_NUM_HC}_hk.csv
# The outputs will then be all the viz stuff. 
# the process downstream of this will gather all individual gc dataframes, merge them, as well as using the gctree 
# API to put all the trees into a forest for any visualizations wanting to be done across many of the trees that have already been featurized. 






# TODO describe
## Post-GCtree processing
#- import idmap, add to GCdf
#- import putative nodes, fill idmap info, add to GCdf
#- calculate total mut#
#- calculate DMS scores
