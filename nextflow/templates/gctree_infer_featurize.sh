#!/bin/bash
set -euo pipefail


# This file is the template for running any particular action
# on the inference, organization with final dataset, node features
# on any one gc. The pipeline will manage the 
# the three parameters below giving you a copy of the full dataset,
# and the gc num that this process is currently running (since it will be in parallel)


# ========================================================
# DO NOT EDIT UNLESS EDITING THE MAIN WORKFLOW AND MODULES


# incase we can multithread inference on a single tree
CPUS=!{task.cpus}

# This variable points towards a df containing a single unique
# germinal center cell type from a specific mouse.
GCDF=!{single_mouse_gc_df}
echo $GCDF

# Strip the .csv suffix for clean def
# Should be annotated-PR-<?>-<?>-<mouse_HC>-<GC_num_HC>-<cell_type_HC>
# TODO remove 'annotated' prefix
GC_DEF="$(basename $GCDF .csv)"

# REQUIRED OUTPUTS
# Whether or not these get populated, 
# this script MUST produce these outputs
# Currently, otherwise the workflow will need some work
#touch ${GC_DEF}-gc-featurized-df.csv    # Featurized df, DMS info etc

# Parse the definition
tokens=(${GC_DEF//\-/ })
PR="${tokens[1]}-${tokens[2]}-${tokens[3]}"
MOUSE=${tokens[4]}
NODE=${tokens[5]}
GC_NUM=${tokens[6]}
CELL_TYPE=${tokens[7]}

# All key row specific results from this template should end up here
OUTDIR="${tokens[1]}${tokens[2]}.${tokens[3]}-${MOUSE}-${NODE}-${GC_NUM}-${CELL_TYPE}"

# Final Kd's from individual (?) mutations
DMS_VSCORES=!{params.dms_vscores}

# Gives you the DMS wild type sites
DMS_SITES=!{params.dms_sites}


# ADD FEATURIZED SEQS TO
mkdir $OUTDIR && cp observed_seqs.csv $OUTDIR

# parameters for hdag mutation models (?)
SUB=!{params.reads_prefix}/!{params.hdag_sub}
MUT=!{params.reads_prefix}/!{params.hdag_mut}
IGK_IDX=!{params.igk_idx}

ls $SUB
ls $MUT

export MPLBACKEND=Agg

# ========================================================
# EDIT BELOW HERE USING PROVIDED VARIABLES ABOVE


# ADD DMS DATA TO OBSERVED BCR SEQS
gcreplay-tools.py featurize-seqs \
    $GCDF \
    --variant_scores ${DMS_VSCORES} \
    --naive_sites ${DMS_SITES} \
    --igk_idx ${IGK_IDX} \
    --output observed_seqs.csv


# check that the file has more than 10 lines (?)
# TODO Should we? I think it fails if not
# if [[ $(wc -l <$GCDF) -ge 15 ]]

# Check that a file is a cell type = GC
# If so we'll be doing the gctree inference.
if [[ "$CELL_TYPE" == "GC" ]];
then

# will not need to do this in the pipeline 

# convert some columns into headers and sequences for a basic fasta format conversion
gcreplay-tools.py gc-df-to-fasta \
    --gc-df $GCDF \
    -o $GC_DEF.fasta
echo \(LOG\) done: created fasta


# make isotype map for the specific gc
gcreplay-tools.py get-columns \
    -df $GCDF \
    -c ID_HK -c isotype_HC \
    -o $GC_DEF.isotypemap
echo \(LOG\) done: extracted isotypes

## First we deduplicate the sequences and convert to phylip alignment format, 
## and also determine sequence abundances. The deduplicate command writes 
## the phylip alignment of unique sequences to stdout (which we redirect to a file here). 
deduplicate $GC_DEF.fasta \
    --root naive \
    --idmapfile $GC_DEF.idmap \
    --abundance_file abundances.csv > deduplicated.phylip
echo \(LOG\) done: deduplicate

# TODO describe
mkconfig deduplicated.phylip dnapars > dnapars.cfg
echo \(LOG\) done: mkconfig

# TODO describe
rm -rf outtree outfile
dnapars < dnapars.cfg > dnapars.log
echo \(LOG\) done: dnapars

cp outfile $OUTDIR
cp abundances.csv $OUTDIR
cp $GC_DEF.isotypemap $OUTDIR
cp $GC_DEF.idmap $OUTDIR

# Run inference on sequences, 
TMPDIR="/tmp"
#mkdir -p ${GC_DEF}-gctree-infer-output/  # all trees
xvfb-run -a gctree infer outfile abundances.csv \
    --idmapfile $GC_DEF.idmap \
    --isotype_mapfile $GC_DEF.isotypemap \
    --chain_split $IGK_IDX \
    --ranking_coeffs 0.1 0.0001 0 \
    --summarize_forest \
    --mutability $MUT \
    --substitution $SUB \
    --root naive \
    --verbose \
    --outbase ${OUTDIR}/gctree \
    | tee gctree.inference.log
echo \(LOG\) done: gctree


# Run will's featurize code
#mkdir -p ${GC_DEF}-featurize-output/     # featurized rank 1 trees? could cobine with below 
xvfb-run -a gcreplay-tools.py featurize-nodes \
    ${OUTDIR}/gctree.inference.1.p \
    ${GC_DEF}.idmap \
    ${DMS_VSCORES} \
    ${DMS_SITES} \
    --igk_idx ${IGK_IDX} \
    --output_dir ${OUTDIR}

mv $GC_DEF.idmap ${OUTDIR}/
echo \(LOG\) done: Viz


# PB/MB Cells
# Just nead to add DMS info to the df.
# else




fi

