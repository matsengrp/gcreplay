#!/bin/bash
set -euo pipefail

CHAIN_LENGTH=50000000
LOG_EVERY=10000

FILENAME=!{observed_seqs_with_time}
BEAST_TEMPLATE=!{beast_template}
OUTDIR=btt-$(basename $FILENAME _with_time.fasta)
mkdir $OUTDIR

beastgen \
    -date_prefix '@' \
    -date_order -1 \
    -D "chain_length=$CHAIN_LENGTH,log_every=$LOG_EVERY" \
    $BEAST_TEMPLATE \
    $FILENAME \
    $OUTDIR/beastgen.xml

beast -working ${OUTDIR}/beastgen.xml
