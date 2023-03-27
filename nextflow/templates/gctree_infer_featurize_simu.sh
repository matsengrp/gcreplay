#!/bin/bash
set -euo pipefail


# =========================================================
# DO NOT* EDIT UNLESS EDITING THE MAIN WORKFLOW AND MODULES

# incase we can multithread inference on a single tree
CPUS=!{task.cpus}

# This is a full path to a simulated event
# Example: taraki-gct-plain-sel-v6-selection-simu-event-1-simu.fasta
SGC_FASTA=!{simulated_gc}
echo $SGC_FASTA

OUTDIR=gct-$(basename $SGC_FASTA .fasta)
mkdir $OUTDIR
cp $SGC_FASTA $OUTDIR/simu.fasta

# Final Kd's from individual (?) mutations
DMS_VSCORES=!{params.dms_vscores}

# Gives you the DMS wild type sites
DMS_SITES=!{params.dms_sites}

# torchdms models
DMS_MULTI_SCORES=!{params.reads_prefix}/!{params.dms_mvscores}
TDMS_MODEL=!{params.reads_prefix}/!{params.tdms_model}
TDMS_MODEL_LINEAR=!{params.reads_prefix}/!{params.tdms_model_lin}

# parameters for hdag mutation models (?)
SUB=!{params.reads_prefix}/!{params.hdag_sub}
MUT=!{params.reads_prefix}/!{params.hdag_mut}
IGK_IDX=!{params.igk_idx}

export MPLBACKEND=Agg
export PYTHONHASHSEED=0

## =========================================================
## EDIT BELOW HERE USING PROVIDED VARIABLES ABOVE FOR GCTREE

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
xvfb-run -a gcreplay-tools.py featurize-nodes \
    ${OUTDIR}/gctree.inference.1.p \
    sim-gc.idmap \
    ${DMS_VSCORES} \
    --multi_variant_scores ${DMS_MULTI_SCORES} \
    --tdms_model ${TDMS_MODEL} \
    --tdms_model_linear ${TDMS_MODEL_LINEAR} \
    --tau 1.0 \
    --tau0 1.0 \
    ${DMS_SITES} \
    --igk_idx ${IGK_IDX} \
    --output_dir ${OUTDIR}

echo \(LOG\) done: Viz
