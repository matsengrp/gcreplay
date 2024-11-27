#!/bin/bash
# -----------------------------------------------------------------------------
# This script is run from the root of the nextflow directory in the
# matsengrp/gcreplay repository, on a gizmo node 
# via `grabnode` command specifying
# 16 cpus and 64 gigs of memory from the Fred Hutch
# computing cluster. 

# nextflow.config specifies the profile being used to
# farm computing steps to the cluster.

# of note, it seems important to specify a the maximum number of forks?
# though I don't have a great understanding of how the nextflow scheduler
# -----------------------------------------------------------------------------


set -euo pipefail

# Load the Nextflow module (if running on rhino/gizmo)
ml nextflow

# Load the Apptainer module (if running on rhino/gizmo with Apptainer)
ml Apptainer

# define the path to the Apptainer binary
export PATH=$APPTAINERROOT/bin/:$PATH
export APPTAINER_TMPDIR=$TMPDIR

# define where you would like the results to be stored
results_dir="results/archive/$(date -I)-full"

# run the Nextflow pipeline
nextflow -log $results_dir/log.log \
        run main.nf \
        --manifest "manifest.csv" \
        --reads_prefix "data/input/" \
        --results $results_dir \
        -profile fred_hutch_rhino \
        -with-trace $results_dir/trace.txt \
        -with-report $results_dir/report.html \
        -resume
