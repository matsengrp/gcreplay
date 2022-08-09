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

# This is a full path to a simulated event
# Example: taraki-gct-plain-sel-v6-selection-simu-event-1-simu.fasta
SGC_FASTA=!{simulated_gc}
echo $SGC_FASTA

#if [[ !{params.test} ]];
#then
#head $SGC_FASTA > grep -v "mrca" > "test-!{simulated_gc}"
#SGC_FASTA="test-!{simulated_gc}"
#fi


OUTDIR=gct-$(basename $SGC_FASTA)
mkdir $OUTDIR
cp $SGC_FASTA $OUTDIR/simu.fasta

# Final Kd's from individual (?) mutations
DMS_VSCORES=!{params.dms_vscores}

# Gives you the DMS wild type sites
DMS_SITES=!{params.dms_sites}


DMS_MULTI_SCORES=!{params.reads_prefix}/!{params.dms_mvscores}
TDMS_MODEL=!{params.reads_prefix}/!{params.tdms_model}
TDMS_MODEL_LINEAR=!{params.reads_prefix}/!{params.tdms_model_lin}

# parameters for hdag mutation models (?)
SUB=!{params.reads_prefix}/!{params.hdag_sub}
MUT=!{params.reads_prefix}/!{params.hdag_mut}
IGK_IDX=!{params.igk_idx}

#ls $SUB
#ls $MUT
#
export MPLBACKEND=Agg
#
## ========================================================
## EDIT BELOW HERE USING PROVIDED VARIABLES ABOVE

## First we deduplicate the sequences and convert to phylip alignment format, 
## and also determine sequence abundances. The deduplicate command writes 
## the phylip alignment of unique sequences to stdout (which we redirect to a file here). 
deduplicate $SGC_FASTA \
    --root naive \
    --idmapfile sim-gc.idmap \
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
cp sim-gc.idmap $OUTDIR

# Run inference on sequences, 
TMPDIR="/tmp"
#mkdir -p ${GC_DEF}-gctree-infer-output/  # all trees
xvfb-run -a gctree infer outfile abundances.csv \
    --idmapfile sim-gc.idmap \
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

cp gctree.inference.log ${OUTDIR}


# Run will's featurize code
#mkdir -p ${GC_DEF}-featurize-output/     # featurized rank 1 trees? could cobine with below 
xvfb-run -a gcreplay-tools.py featurize-nodes \
    ${OUTDIR}/gctree.inference.1.p \
    ${GC_DEF}.idmap \
    ${DMS_VSCORES} \
    --multi_variant_scores ${DMS_MULTI_SCORES} \
    --tdms_model ${TDMS_MODEL} \
    --tdms_model_linear ${TDMS_MODEL_LINEAR} \
    --tau 1.0 \
    --tau0 1.0 \
    ${DMS_SITES} \
    --igk_idx ${IGK_IDX} \
    --output_dir ${OUTDIR}

mv $GC_DEF.idmap ${OUTDIR}/
echo \(LOG\) done: Viz
#
#
## PB/MB Cells
## Just nead to add DMS info to the df.
## else
#
#
#
#
#fi

