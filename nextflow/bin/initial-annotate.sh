#!/bin/bash
# example usage:
#  ./datascripts/meta/taraki-gctree-2021-10/initial-annotate.sh /fh/fast/matsen_e/dralph/partis/merged_bcrs.fasta $fs/partis/tmp/gcreplay-test /home/dralph/work/partis/datascripts/meta/taraki-gctree-2021-10/germlines

bin=./bin/partis

ifn=$1
odir=$2
gldir=$3

for rstr in engrd all; do
    common="--paired-loci --no-pairing-info --species mouse --airr-output"
    if [ "$rstr" == "engrd" ]; then
	common="$common --no-insertions-or-deletions --leave-default-germline --initial-germline-dir $gldir"
    fi
    for act in cache-parameters annotate; do
	mkdir -p $odir/$rstr
	cmd="$bin $act $common --infname $ifn --paired-outdir $odir/$rstr"
	echo $cmd
	$cmd >$odir/$rstr/$act.log
    done
done
