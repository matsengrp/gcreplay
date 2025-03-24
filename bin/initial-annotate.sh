#!/bin/bash
# example usage:
#  ./datascripts/meta/taraki-gctree-2021-10/initial-annotate.sh /fh/fast/matsen_e/dralph/partis/merged_bcrs.fasta $fs/partis/tmp/gcreplay-test /home/dralph/work/partis/datascripts/meta/taraki-gctree-2021-10/germlines

bin=/partis/bin/partis

INPUT_FASTA=$1
OUTPUT_DIR=$2
GERMLINE_DIR=$3

for rstr in engrd; do
   
   common="--random-seed 23 --paired-loci --no-pairing-info --species mouse --airr-output --extra-annotation-columns seqs_aa:n_mutations"

   # if we're looking only at the engineered sequences, use the flags
   if [ "$rstr" == "engrd" ]; then
	    common="$common --no-insertions-or-deletions --leave-default-germline --initial-germline-dir $GERMLINE_DIR"
   fi

   for act in cache-parameters annotate; do
	    mkdir -p $OUTPUT_DIR/$rstr
	    cmd="$bin $act $common --infname $INPUT_FASTA --paired-outdir $OUTPUT_DIR/$rstr"
	    echo $cmd
	    $cmd >$OUTPUT_DIR/$rstr/$act.log
   done
done