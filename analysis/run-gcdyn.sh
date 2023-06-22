#!/bin/bash

# TEST ./projects/cf-gcdyn.py --actions simu:dl-infer --n-replicates 1 --label test --version v0 --dry --test --n-trials 6 --n-sub-procs 3

bin=./projects/cf-gcdyn.py
common="--actions simu --dry --n-max-procs 5 --base-outdir  /fh/fast/matsen_e/dralph/partis/gcdyn" #/fh/local/dralph/partis/gcdyn"
# echo ./projects/cf-gcdyn.py --actions simu --version v1 --birth-response-list constant:soft-relu:soft-relu:soft-relu --xscale-list 1:1:2:3 --zip-vars birth-response:xscale $commone
xscales=0.1:1:5:10; xshifts=-5:-2:-1:0:1:2:5
# common="--actions simu --final-plot-xvar xshift --dry --n-replicates 2"
# echo $bin $common --version v2 --birth-response-list soft-relu --xscale-list $xscales --xshift-list=$xshifts
# echo $bin $common --version v3-carry-cap --birth-response-list sigmoid --xscale-list $xscales --xshift-list=$xshifts
# echo $bin $common --version v4-carry-cap --birth-response-list sigmoid --xscale-list 0.5:0.75:1:1.25:1.5 --xshift-list=0:1:2:5:8:15
# common="--actions simu:merge-simu:process --final-plot-xvar xshift --dry"
# echo $bin $common --version v6 --birth-response-list sigmoid --xscale-list 0.75:1:1.25 --xshift-list=2:5 --n-replicates 20
# echo $bin $common --label vs-n-trees --version v0 --birth-response-list sigmoid --carry-cap-list 150 --xscale-list 0.5,1.5 --xshift-list 2 --n-replicates 2 --n-trials-list 100:1000 --simu-extra-args=\"--n-sub-procs 10\" --perf-metrics xscale --plot-metrics dl-infer
# echo $bin $common --label vs-n-trees --version v1 --birth-response-list sigmoid --carry-cap-list 150 --xscale-list 1 --xshift-list 2 --n-replicates 2 --n-trials-list 100:500 --simu-extra-args=\"--n-sub-procs 10\" --perf-metrics xscale --plot-metrics dl-infer
# echo $bin $common --label vs-n-leaves --version v2 --carry-cap-list 150 --xscale-list 0.5,0.75,1,1.25,1.5 --xshift-list 2 --n-replicates 2 --n-seqs-list 50:75:100:125 --n-trials-list 100:1000 --simu-extra-args=\"--n-sub-procs 10\" --perf-metrics all-dl --model-size-list small:tiny --plot-metrics dl-infer --final-plot-xvar n-seqs --pvks-to-plot 1000
echo $bin $common --label vs-n-trees --version v4 --carry-cap-list 150 --xscale-list 0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5 --xshift-list 2 --n-replicates 3 --n-trials-list 500:5000:50000 --test-xscale-values-list 0.7,0.8,0.9,1,1.1,1.2,1.3 --simu-extra-args=\"--n-sub-procs 10\" --perf-metrics all-dl --plot-metrics dl-infer --final-plot-xvar n-trials
echo "--perf-metrics xscale-test-variance --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5:10:50:100"

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

