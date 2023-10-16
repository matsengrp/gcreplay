#!/bin/bash

# TEST ./projects/cf-gcdyn.py --actions simu:dl-infer --n-replicates 1 --label test --version v0 --dry --test --n-trials 6 --n-sub-procs 3

bin=./projects/cf-gcdyn.py
common="--actions simu --dry --n-max-procs 5 --n-sub-procs 100 --base-outdir  /fh/fast/matsen_e/dralph/partis/gcdyn" #/fh/local/dralph/partis/gcdyn"
# echo ./projects/cf-gcdyn.py --actions simu --version v1 --birth-response-list constant:soft-relu:soft-relu:soft-relu --xscale-list 1:1:2:3 --zip-vars birth-response:xscale $commone
xscales=0.1:1:5:10; xshifts=-5:-2:-1:0:1:2:5
# common="--actions simu --final-plot-xvar xshift --dry --n-replicates 2"
# echo $bin $common --version v2 --birth-response-list soft-relu --xscale-list $xscales --xshift-list=$xshifts
# echo $bin $common --version v3-carry-cap --birth-response-list sigmoid --xscale-list $xscales --xshift-list=$xshifts
# echo $bin $common --version v4-carry-cap --birth-response-list sigmoid --xscale-list 0.5:0.75:1:1.25:1.5 --xshift-list=0:1:2:5:8:15
# common="--actions simu:merge-simu:process --final-plot-xvar xshift --dry"
# echo $bin $common --version v6 --birth-response-list sigmoid --xscale-list 0.75:1:1.25 --xshift-list=2:5 --n-replicates 20
# echo $bin $common --label vs-n-trees --version v0 --birth-response-list sigmoid --carry-cap-list 150 --xscale-list 0.5,1.5 --xshift-list 2 --n-replicates 2 --n-trials-list 100:1000 --perf-metrics xscale --plot-metrics dl-infer
# echo $bin $common --label vs-n-trees --version v1 --birth-response-list sigmoid --carry-cap-list 150 --xscale-list 1 --xshift-list 2 --n-replicates 2 --n-trials-list 100:500 --perf-metrics xscale --plot-metrics dl-infer
# echo $bin $common --label vs-n-leaves --version v2 --carry-cap-list 150 --xscale-list 0.5,0.75,1,1.25,1.5 --xshift-list 2 --n-replicates 2 --n-seqs-list 50:75:100:125 --n-trials-list 100:1000 --perf-metrics all-dl --model-size-list small:tiny --plot-metrics dl-infer --final-plot-xvar n-seqs --pvks-to-plot 1000
# echo $bin $common --label vs-n-trees --version v6 --carry-cap-list 150 --xscale-list 0.5,1,1.5,2,2.5,3,3.5,4,4.5,5 --xshift-list=-0.5,0,0.5,1,1.5,2,2.5,3 --n-replicates 2 --n-trials-list 500:5000:50000:500000 --test-xscale-values-list 1.5,2,2.5,3,3.5,4 --test-xshift-values-list 0.5,1,1.5,2,2.5 --simu-extra-args=\"--n-max-procs 20\" --perf-metrics all-test-dl --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5:10:50:100
# echo $bin $common --label vs-n-trees --version v7 --carry-cap-list 150 --xscale-list 2 --xshift-list=-0.5,3 --n-replicates 2 --n-trials-list 500:5000:50000 --simu-extra-args=\"--n-max-procs 20\" --perf-metrics all-test-dl --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5:10:50:100
# echo ./projects/cf-gcdyn.py --actions dl-infer --n-max-procs 2 --base-outdir /fh/fast/matsen_e/dralph/partis/gcdyn --label test-simu-test --version v1 --carry-cap-list 150 --xscale-values-list 2 --xshift-range-list=-0.5,3 --n-replicates 2 --n-trials-list 500:5000 --simu-extra-args="--n-max-procs 20 --n-trees-per-param-set 10" --dl-extra-args="--no-shuffle" --perf-metrics all-test-dl --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5 --params-to-predict xshift --dry
# echo ./projects/cf-gcdyn.py --actions dl-infer --n-max-procs 4 --base-outdir /fh/fast/matsen_e/dralph/partis/gcdyn --label vary-trees-per-param-set --version v0 --carry-cap-list 150 --xscale-values-list 2 --xshift-range-list=-0.5,3 --n-replicates 2 --n-trials-list 500:5000 --n-trees-per-param-set-list 10:50 --simu-extra-args="--n-max-procs 20" --dl-extra-args="--no-shuffle" --perf-metrics all-test-dl --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5 --params-to-predict xshift --dry
echo $bin $common --label bundle-xshift --version v8 --carry-cap-list 150 --xscale-values-list 2 --xshift-range-list=-0.5,3 --n-replicates 1 --n-trials-list 5000:50000 --simu-bundle-size 50 --simu-extra-args=\"--n-max-procs 20\" --perf-metrics all-test-dl --plot-metrics group-expts --final-plot-xvar n-trees-per-expt --n-trees-per-expt-list 1:5 --params-to-predict xshift
# NOTE made bundle size 250 for 500k sample (could use zip vars for this, but also want to zip n-trials with epochs)
echo $bin $common --label bundle-xscale-xshift --version v0 --carry-cap-list 150 --xscale-range-list 0.5,5 --xshift-range-list=-0.5,3 --n-replicates 1 --n-trials-list 5000:50000:500000 --epochs-list 1000:5000:1000 --zip-vars n-trials:epochs --simu-bundle-size 50 --simu-extra-args=\"--n-max-procs 20\" --params-to-predict xscale:xshift --dropout-rate-list 0:0.2:0.5 --learning-rate-list 0.001:0.01 --ema-momentum-list 0.9:0.99
echo $bin $common --label new-constraints --version v0 --carry-cap-list 150 --xscale-range-list 0.5,5 --xshift-range-list=-0.5,3 --yscale-range-list 1,50 --initial-birth-rate-range-list 2,10 --n-replicates 1 --n-trials-list 100 --time-to-sampling-values-list 5 --simu-extra-args=\"--n-max-procs 20\" --params-to-predict xscale:xshift:yscale
# --params-to-predict xshift
# xshift 1  update:maybe not
# include larger xscale values like 5-10
# also try to infer xshift

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

