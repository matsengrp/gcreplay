#!/bin/bash
set -euo pipefail

# CHAIN_LENGTH=50000000
# LOG_EVERY=10000
CHAIN_LENGTH=5000
LOG_EVERY=100

CPUS=!{task.cpus}
FILENAME=!{observed_seqs_with_time}
BEAST_TEMPLATE=!{beast_template}

# _with_time is added before, strip that
BASENAME=$(basename $FILENAME _with_time.fasta)

# remove the 'annotated' prefix from the fasta if it exists.
OUTDIR=btt-${BASENAME#"annotated-"}

# all IO from beast will go here for downstream conversion
mkdir $OUTDIR

beastgen \
    -date_prefix '@' \
    -date_order -1 \
    -D "chain_length=$CHAIN_LENGTH,log_every=$LOG_EVERY" \
    $BEAST_TEMPLATE \
    $FILENAME \
    $OUTDIR/beastgen.xml

beast -threads ${CPUS} -prefix beast -working $OUTDIR/beastgen.xml
