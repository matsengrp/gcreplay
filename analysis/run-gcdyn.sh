#!/bin/bash

vsn=test-v0
ddir=/fh/fast/matsen_e/dralph/gcdyn/gcreplay-observed
odir=/fh/fast/matsen_e/dralph/gcdyn/$vsn
gcrdir=projects/gcreplay/analysis
gcddir=projects/gcdyn

echo ./$gcrdir/gcdyn-plot.py --data-dir $ddir --simu-dir $odir/simu --outdir $odir/plots
# echo conda activate gcdyn
# echo python $gcddir/scripts/multi-simulation.py --outdir $odir/simu --n-seqs 70 --n-trials 1
