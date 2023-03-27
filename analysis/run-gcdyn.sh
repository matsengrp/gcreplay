#!/bin/bash

echo ./projects/cf-gcdyn.py --actions simu --version v1 --birth-response-list constant:soft-relu:soft-relu:soft-relu --xscale-list 1:1:2:3 --zip-vars birth-response:xscale --n-replicates 2
echo ./projects/cf-gcdyn.py --actions simu:process --version v2 --birth-response-list soft-relu --xscale-list 0.1:1:10 --xshift-list=-5:-1:0:1:5 --n-replicates 2 --dry
echo ./projects/cf-gcdyn.py --actions simu:process --version v3 --birth-response-list sigmoid --xscale-list 0.1:1:10 --xshift-list=-5:-1:0:1:5 --n-replicates 2 --dry
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

