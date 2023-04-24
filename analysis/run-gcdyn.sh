#!/bin/bash

bin=./projects/cf-gcdyn.py
common="--actions simu:process --final-plot-xvar xshift --dry --n-replicates 2"
# echo ./projects/cf-gcdyn.py --actions simu --version v1 --birth-response-list constant:soft-relu:soft-relu:soft-relu --xscale-list 1:1:2:3 --zip-vars birth-response:xscale $commone
xscales=0.1:1:5:10; xshifts=-5:-2:-1:0:1:2:5
# echo $bin $common --version v2 --birth-response-list soft-relu --xscale-list $xscales --xshift-list=$xshifts
echo $bin $common --version v3 --birth-response-list sigmoid --xscale-list $xscales --xshift-list=$xshifts
echo $bin $common --version v4 --birth-response-list sigmoid --xscale-list 0.5:0.75:1:1.25:1.5 --xshift-list=0:1:2:5:8:15
exit 0

vsn=test-v2
ddir=/fh/fast/matsen_e/dralph/gcdyn/gcreplay-observed
odir=/fh/fast/matsen_e/dralph/gcdyn/$vsn
gcrdir=projects/gcreplay/analysis
gcddir=projects/gcdyn

echo conda activate gcdyn
echo python $gcddir/scripts/multi-simulation.py --outdir $odir/simu --n-seqs 70 --n-trials 51
echo ./$gcrdir/gcdyn-plot.py --data-dir $ddir --simu-dir $odir/simu --outdir $odir

# simulation scan
#  - move into scanplot (i.e. make a new script projecrs/cf-xxx.py
#  - keep track of # of tries and time per family

