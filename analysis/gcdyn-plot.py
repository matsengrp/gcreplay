#!/usr/bin/env python2
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import glob
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = '%s/work/partis' % os.getenv('HOME') #os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import plotting
import lbplotting
import hutils

# ----------------------------------------------------------------------------------------
colors = {'data' : plotting.default_colors[0],
          'simu' : plotting.default_colors[1]}

# ----------------------------------------------------------------------------------------
def parse_name(fn):  # convert .fa name to mouse #, etc
    fn = os.path.basename(fn)
    # fn = 'annotated-PR-1-10-11-LA-122-GC.fasta'
    nlist = fn.split('-')
    if len(nlist) != 8:
        raise Exception('fn should have 7 - chars, but got %d: %s' % (fn.count('-'), fn))
    assert nlist[:2] == ['annotated', 'PR']
    flowcell = '-'.join(nlist[2:4])
    mouse = int(nlist[4])
    ln_loc = nlist[5]
    ln_id = nlist[6]
    assert nlist[7] == 'GC.fasta'
    return mouse, flowcell, ln_loc, ln_id

# ----------------------------------------------------------------------------------------
def bstr(mouse, flowcell, ln_loc, ln_id):
    return '-'.join([str(mouse), flowcell, ln_loc, ln_id])
# ----------------------------------------------------------------------------------------
def filter_mice(in_fns):
    skipped_mice, kept_mice = set(), set()
    fn_fns = []
    for fasta_path in in_fns:
        mouse, flowcell, ln_loc, ln_id = parse_name(fasta_path)
        if args.mice is not None and mouse not in args.mice:
            skipped_mice.add(mouse)
            continue
        kept_mice.add(mouse)
        fn_fns.append(fasta_path)
    if len(skipped_mice) > 0:
        print '    skipped %d mice: %s' % (len(skipped_mice), ' '.join(str(s) for s in sorted(skipped_mice)))
    print '    kept %d samples from %d mice: %s' % (len(fn_fns), len(kept_mice), ' '.join(str(s) for s in sorted(kept_mice)))
    return fn_fns

# ----------------------------------------------------------------------------------------
def abfn(tlab):
    return '%s/%s/abundances.csv' % (args.outdir, tlab)

# ----------------------------------------------------------------------------------------
def calc_abdn(indir, label, is_simu):
    fn_fns = glob.glob('%s/*.fasta'%indir)
    if not is_simu:
        fn_fns = filter_mice(fn_fns)
    if len(fn_fns) == 0:
        raise Exception('no fasta files in dir %s' % indir)
    cmd = 'python %s/scripts/abundance.py %s --min-seqs 70 --max-seqs 70 --outdir %s' % (args.gcdyn_dir, ' '.join(fn_fns), os.path.dirname(abfn(label)))
    utils.simplerun(cmd)

# ----------------------------------------------------------------------------------------
def hargs(htmp):
    xbounds = [0.5, htmp.xmax+0.5]
    xticks = list(range(1, int(htmp.xmax)+1, 2))
    ybounds = [0.9 * htmp.get_minimum(exclude_empty=True), 1.1 * htmp.get_maximum()]
    yticks, yticklabels = plotting.get_auto_y_ticks(ybounds[0], ybounds[1])
    ybounds = [yticks[0], yticks[-1]]
    return xbounds, ybounds, xticks, yticks, yticklabels

# ----------------------------------------------------------------------------------------
def plot(plotdir, label):
    with open(abfn(label)) as afile:
        reader = csv.DictReader(afile)
        plotvals = {k : {} for k in reader.fieldnames if k!=''}
        for line in reader:
            abn = int(line[''])  # i can't figure out how to set this column label in the other script
            for bn in plotvals:  # <bn> is GC name
                plotvals[bn][abn] = int(line[bn])

    utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

    max_vals = []
    distr_hists = []
    for bn in plotvals:
        max_vals.append(max([a for a, n in plotvals[bn].items() if n>0]))
        htmp = hutils.make_hist_from_dict_of_counts(plotvals[bn], 'int', bn)
        if args.max_gc_plots is not None and len(max_vals) <= args.max_gc_plots:
            fn = htmp.fullplot(plotdir, 'abdn-distr-gc-%s'%bn, pargs={'square_bins' : True, 'errors' : False, 'color' : colors[label]}, fargs={'title' : bn, 'xlabel' : 'abundance', 'ylabel' : 'counts'})
            lbplotting.add_fn(fnames, fn)
        distr_hists.append(htmp)

    hmax = hutils.make_hist_from_list_of_values(max_vals, 'int', 'max-abdn')
    xbounds, ybounds, xticks, yticks, yticklabels = hargs(hmax)
    hmax.title = label
    hmax.xtitle = 'max abundance in GC'
    fn = hmax.fullplot(plotdir, 'max-abdn', pargs={'square_bins' : True, 'remove_empty_bins' : True, 'errors' : False, 'color' : colors[label]},
                       fargs={'xbounds' : xbounds, 'ybounds' : ybounds, 'yticks' : yticks, 'yticklabels' : yticklabels,
                              'xlabel' : hmax.xtitle, 'ylabel' : 'counts', 'log' : 'y', 'xticks' : xticks})
    # fnames[0].append(fn)

    mean_hdistr = plotting.make_mean_hist(distr_hists)
    mean_hdistr.title = '%s (%d GCs)' % (label, len(distr_hists))
    mean_hdistr.xtitle = 'abundance'
    mean_hdistr.ytitle = 'N seqs\nmean+/-std, %d GCs' % len(distr_hists)
    xbounds, ybounds, xticks, yticks, yticklabels = hargs(mean_hdistr)
    fn = mean_hdistr.fullplot(plotdir, 'abdn', pargs={'remove_empty_bins' : True, 'color' : colors[label]},
                              fargs={'xbounds' : xbounds, 'ybounds' : ybounds, 'yticks' : yticks, 'yticklabels' : yticklabels,
                                     'xlabel' : mean_hdistr.xtitle, 'ylabel' : mean_hdistr.ytitle, 'log' : 'y', 'xticks' : xticks})
    # fnames[0].append(fn)
    return {'distr' : mean_hdistr, 'max' : hmax}

# ----------------------------------------------------------------------------------------
def compare_plots(hname, plotdir, hists, labels):
    xbounds, ybounds, xticks, yticks, yticklabels = hargs(hists[0])
    ytitle = hists[0].ytitle
    if hname == 'distr':
        ytitle = '%s\nmean +/- std' % hists[0].ytitle.split('\n')[0]
    fn = plotting.draw_no_root(None, plotdir=plotdir, plotname='%s-abdn'%hname, more_hists=hists, log='y', xtitle=hists[0].xtitle, ytitle=ytitle,
                               bounds=xbounds, ybounds=ybounds, xticks=xticks, yticks=yticks, yticklabels=yticklabels, errors=hname!='max', square_bins=hname=='max', linewidths=[4, 3], plottitle='',
                               alphas=[0.6, 0.6], colors=[colors[l] for l in labels])
    fnames[0].append(fn)

# ----------------------------------------------------------------------------------------
ustr = """
./scripts/plot-abdn.py --data-dir <path with fastas> --simu-dir <path with fastas>
"""
parser = argparse.ArgumentParser(usage=ustr)
parser.add_argument('--data-dir')
parser.add_argument('--simu-dir')
parser.add_argument('--outdir')
parser.add_argument('--mice', default=[1, 2, 3, 4, 5, 6], help='restrict to these mouse numbers')
parser.add_argument('--is-simu', action='store_true')
parser.add_argument('--gcdyn-dir', default='%s/work/partis/projects/gcdyn'%os.getenv('HOME'))
parser.add_argument('--max-gc-plots', type=int, default=2, help='only plot individual (per-GC) plots for this  many GCs')
args = parser.parse_args()

dlabels = []
if args.data_dir is not None:
    dlabels.append([args.data_dir, 'data'])
if args.simu_dir is not None:
    dlabels.append([args.simu_dir, 'simu'])

hclists, fnames = {'distr' : [], 'max' : []}, [[], [], []]
for idir, tlab in dlabels:
    calc_abdn(idir, tlab, is_simu=tlab=='simu')
    lhists = plot('%s/plots/%s'%(args.outdir, tlab), tlab)
    for hn in lhists:
        hclists[hn].append(lhists[hn])

cfpdir = '%s/plots/comparisons' % args.outdir
utils.prep_dir(cfpdir, wildlings=['*.csv', '*.svg'])
for hname, hlist in hclists.items():
    compare_plots(hname, cfpdir, hlist, [l for _, l in dlabels])
plotting.make_html(args.outdir+'/plots', fnames=fnames)
