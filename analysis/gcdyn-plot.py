#!/usr/bin/env python2
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import glob
import argparse
import colored_traceback.always
import json

# if you move this script, you'll need to change this method of getting the imports
partis_dir = '%s/work/partis' % os.getenv('HOME') #os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import plotting
import lbplotting
import hutils
from hist import Hist

# ----------------------------------------------------------------------------------------
colors = {'data' : plotting.default_colors[1],
          'simu' : plotting.default_colors[0]}
pltlabels = {'hdists' : 'root-tip dist', 'max-abdn-shm' : 'median SHM of seqs\nw/max abundance'}

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
def abfn(tlab, abtype='abundances'):
    return '%s/%s/%s.csv' % (args.outdir, tlab, abtype)

# ----------------------------------------------------------------------------------------
def parse_fastas(indir, label, is_simu):
    fn_fns = glob.glob('%s/%s*.fasta' % (indir, 'seqs_' if is_simu else 'annotated-'))
    if not is_simu:
        fn_fns = filter_mice(fn_fns)
    if len(fn_fns) == 0:
        raise Exception('no fasta files in dir %s' % indir)
    cmd = 'python %s/scripts/abundance.py %s --min-seqs 70 --max-seqs 70 --outdir %s' % (args.gcdyn_dir, ' '.join(fn_fns), os.path.dirname(abfn(label)))
    utils.simplerun(cmd)

# ----------------------------------------------------------------------------------------
def hargs(hlist):
    xmax = max(h.xmax for h in hlist)
    xbounds = [0.5, xmax+0.5]
    xticks = list(range(1, int(xmax)+1, 2 if xmax<10 else 5))
    ymin = min(h.get_minimum(exclude_empty=True) for h in hlist)
    ymax = max(h.get_maximum() for h in hlist)
    ybounds = [0.9 * ymin, 1.1 * ymax]
    yticks, yticklabels = plotting.get_auto_y_ticks(ybounds[0], ybounds[1])
    ybounds = [yticks[0], yticks[-1]]
    return xbounds, ybounds, xticks, yticks, yticklabels

# ----------------------------------------------------------------------------------------
def plot(plotdir, label, abtype):
    # read values from csv file
    with open(abfn(label, abtype=abtype)) as afile:
        reader = csv.DictReader(afile)
        if abtype == 'abundances':
            plotvals = {k : {} for k in reader.fieldnames if k!=''}
            for line in reader:
                abn = int(line[''])  # i can't figure out how to set this column label in the other script
                for bn in plotvals:  # <bn> is GC name
                    plotvals[bn][abn] = int(line[bn])
        else:
            plotvals = {}
            for line in reader:
                plotvals[line['fbase']] = [int(s) for s in line['vlist'].split(':')]

    # collect values from each GC
    max_vals = []
    distr_hists = []
    for bn in plotvals:
        if abtype == 'abundances':
            max_vals.append(max([a for a, n in plotvals[bn].items() if n>0]))
            htmp = hutils.make_hist_from_dict_of_counts(plotvals[bn], 'int', bn)
        else:
            htmp = hutils.make_hist_from_list_of_values(plotvals[bn], 'int', bn) #, xmin_force=-0.5 if abtype=='max-abdn-shm' else None)
        if len(distr_hists) < args.max_gc_plots:
            nstxt = '%d seqs'%htmp.integral(True, multiply_by_bin_center=abtype=='abundances')
            mvtxt = 'mean %.1f' % htmp.get_mean()
            fn = htmp.fullplot(plotdir, '%s-distr-gc-%s'%(abtype, bn), pargs={'square_bins' : True, 'errors' : False, 'color' : colors[label]},
                               fargs={'title' : bn, 'xlabel' : pltlabels.get(abtype, abtype), 'ylabel' : 'counts', 'log' : 'y', 'title' : '%s: %s'%(label, bn)}, texts=[[0.6, 0.8, nstxt], [0.6, 0.75, mvtxt]])
            lbplotting.add_fn(fnames, fn, new_row=len(distr_hists)==0)
        distr_hists.append(htmp)

    hmax = None
    if abtype == 'abundances':
        hmax = hutils.make_hist_from_list_of_values(max_vals, 'int', 'max-abdn')
        xbounds, ybounds, xticks, yticks, yticklabels = hargs([hmax])
        hmax.title = '%s (mean %.1f)' % (label, hmax.get_mean())
        hmax.xtitle = 'max abundance in GC'

    # plot mean distribution over GCs
    mean_hdistr = plotting.make_mean_hist(distr_hists)
    mean_hdistr.title = '%s (mean %.1f)' % (label, mean_hdistr.get_mean())
    mean_hdistr.xtitle = pltlabels.get(abtype, abtype)
    mean_hdistr.ytitle = 'N seqs\nmean+/-std, %d GCs' % len(distr_hists)
    xbounds, ybounds, xticks, yticks, yticklabels = hargs([mean_hdistr])

    return {'distr' : mean_hdistr, 'max' : hmax}

# ----------------------------------------------------------------------------------------
# NOTE may be better to eventually switch to optimal transport rather than this weighted average/center of mass approach
# NOTE may also/instead want to use log of y difference (The max abundance seems ok, but the high tail of the abundance distr is getting totally washed out/overwhelmed by the low end)
def hist_distance(h1, h2, dbgstr='hist', weighted=False, debug=False):
    if debug:
        print '    %s distance%s:' % (dbgstr, ' (weighted)' if weighted else '')
        print '      xval     v1      v2    abs diff'
    xvals = sorted(set(x for h in [h1, h2] for x in h.get_bin_centers()))
    dvals = []
    for xval in xvals:
        ib1, ib2 = [h.find_bin(xval) for h in [h1, h2]]
        if not h1.is_overflow(ib1) and not h2.is_overflow(ib2):
            if h1.low_edges[ib1] != h2.low_edges[ib2]:
                raise Exception('h1 low edge %.3f (ibin %d) doesn\'t equal h2 low edge %.3f (ibin %d)' % (h1.low_edges[ib1], ib1, h2.low_edges[ib2], ib2))
        v1, v2 = [0. if h.is_overflow(i) else h.bin_contents[i] for h, i in zip([h1, h2], [ib1, ib2])]  # NOTE ignores contents of over and underflow bins, which isn't great, but there shouldn't be any?
        dvals.append((xval if weighted else 1) * abs(v1 - v2))
        if debug:
            def fstr(v): return utils.color('blue', '-', width=6) if v==0 else '%6.2f'%v
            print '      %3.0f  %s  %s  %s' % (xval, fstr(v1), fstr(v2), fstr(dvals[-1]))
    if debug:
        print '    %.1f' % sum(dvals)
    return sum(dvals)

# ----------------------------------------------------------------------------------------
def compare_plots(hname, plotdir, hists, labels, abtype, diff_vals):
    xbounds, ybounds, xticks, yticks, yticklabels = hargs(hists) if abtype=='abundances' else (None, None, None, None, None)
    ytitle = hists[0].ytitle
    if hname == 'distr':
        ytitle = '%s\nmean +/- std' % hists[0].ytitle.split('\n')[0]
    fn = plotting.draw_no_root(None, plotdir=plotdir, plotname='%s-%s'%(hname, abtype), more_hists=hists, log='y' if abtype=='abundances' else '', xtitle=hists[0].xtitle, ytitle=ytitle,
                               bounds=xbounds, ybounds=ybounds, xticks=xticks, yticks=yticks, yticklabels=yticklabels, errors=hname!='max', square_bins=hname=='max', linewidths=[4, 3], plottitle='',
                               alphas=[0.6, 0.6], colors=[colors[l] for l in labels], translegend=[-0.2, 0], write_csv=True, hfile_labels=labels)
    fnames[0].append(fn)

    hdict = {l : h for l, h in zip(labels, hists)}
    if abtype != 'max-abdn-shm':  # min/max aren't the same for all hists, so doesn't work yet
        dname = '%s-%s'%(hname, abtype)
        diff_vals[dname] = hist_distance(hdict['simu'], hdict['data'], weighted='abundances' in abtype, dbgstr=dname)

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
parser.add_argument('--max-gc-plots', type=int, default=3, help='only plot individual (per-GC) plots for this  many GCs')
args = parser.parse_args()

dlabels = []
if args.data_dir is not None:
    dlabels.append([args.data_dir, 'data'])
if args.simu_dir is not None:
    dlabels.append([args.simu_dir, 'simu'])

abtypes = ['abundances', 'hdists', 'max-abdn-shm']
hclists = {t : {'distr' : [], 'max' : []} for t in abtypes}
fnames = [[], [], []]
for idir, tlab in dlabels:
    parse_fastas(idir, tlab, is_simu=tlab=='simu')
    utils.prep_dir('%s/plots/%s'%(args.outdir, tlab), wildlings=['*.csv', '*.svg'])
    for abtype in abtypes:
        lhists = plot('%s/plots/%s'%(args.outdir, tlab), tlab, abtype)
        for hn in lhists:
            hclists[abtype][hn].append(lhists[hn])

cfpdir = '%s/plots/comparisons' % args.outdir
utils.prep_dir(cfpdir, wildlings=['*.csv', '*.svg'])
diff_vals = {}
for abtype in abtypes:
    for hname, hlist in hclists[abtype].items():
        if hlist.count(None) == len(hlist):
            continue
        compare_plots(hname, cfpdir, hlist, [l for _, l in dlabels], abtype, diff_vals)

plotting.make_html(args.outdir+'/plots', fnames=fnames)
dfn = '%s/diff-vals.yaml' % args.outdir
with open(dfn, 'w') as dfile:
    json.dump(diff_vals, dfile)
