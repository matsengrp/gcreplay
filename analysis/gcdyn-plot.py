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
pltlabels = {'hdists' : 'root-tip dist', 'max-abdn-shm' : 'SHM in most\nabundant seq'}  # 'median SHM of seqs\nw/max abundance'}

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
    if os.path.exists(abfn(label)):
        print'    %s abundance exists, skipping: %s' % (label, abfn(label))
        return
    if is_simu:
        workdir = '%s/work' % args.outdir
        seqfos = utils.read_fastx('%s/seqs.fasta'%indir)
        naive_seq = None  # collect families from single fasta output file (this is annoying but i think still better than going back to writing all individual fastas for each family)
        family_fos = {}
        for sfo in seqfos:
            if sfo['name'] == 'naive':
                if naive_seq is None:
                    naive_seq = sfo['seq']
                assert sfo['seq'] == naive_seq
                continue
            fid, sid = sfo['name'].split('-')
            if fid not in family_fos:
                family_fos[fid] = []
            family_fos[fid].append(sfo)
        for fid in sorted(family_fos, key=int):
            utils.write_fasta('%s/seqs_%d.fasta'%(workdir, int(fid)), [{'name' : 'naive', 'seq' : naive_seq}] + family_fos[fid])
        fn_fns = glob.glob('%s/seqs_*.fasta' % workdir)
    else:
        fn_fns = glob.glob('%s/annotated-*.fasta' % indir)
    if not is_simu:
        fn_fns = filter_mice(fn_fns)
    if len(fn_fns) == 0:
        raise Exception('no fasta files in dir %s' % indir)
    cmd = 'python %s/scripts/abundance.py %s%s --max-seqs 70 --outdir %s' % (args.gcdyn_dir, ' '.join(fn_fns), '' if is_simu else ' --min-seqs %d'%args.min_seqs_per_gc, os.path.dirname(abfn(label)))
    utils.simplerun(cmd)
    if is_simu:
        utils.clean_files(fn_fns)
        os.rmdir(workdir)

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
def plot_abdn_stuff(plotdir, label, abtype):
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
    for bn in plotvals:  # bn is the base name of the fasta file (i.e. gc name)
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
        hmax.title = '%s (%d trees)' % (label, len(distr_hists))
        hmax.xtitle = 'max abundance in GC'

    # plot mean distribution over GCs
    mean_hdistr = plotting.make_mean_hist(distr_hists)
    mean_hdistr.title = '%s (%d trees)' % (label, len(distr_hists))
    mean_hdistr.xtitle = pltlabels.get(abtype, abtype)
    mean_hdistr.ytitle = 'N seqs in bin\nmean+/-std (over GCs)'
    xbounds, ybounds, xticks, yticks, yticklabels = hargs([mean_hdistr])

    return {'distr' : mean_hdistr, 'max' : hmax}

# ----------------------------------------------------------------------------------------
def get_affinity_plots(label):
    # ----------------------------------------------------------------------------------------
    def get_plotvals(label, dendro_trees, affy_vals):
        def nstr(node):
            return 'leaf' if node.is_leaf() else 'internal'
        plotvals = {k : [] for k in ['leaf', 'internal']}
        n_missing, n_tot = {'internal' : [], 'leaf' : []}, {'internal' : [], 'leaf' : []}
        for dtree in dendro_trees:
            for node in dtree.preorder_node_iter():
                name = node.taxon.label
                n_tot[nstr(node)].append(name)
                if affy_vals.get(name) is None:
                    n_missing[nstr(node)].append(name)
                    continue
                plotvals[nstr(node)].append(affy_vals[name])
        if sum(len(l) for l in n_missing.values()) > 0:
            print '      %s missing/none affinity values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal']))
        return plotvals
    # ----------------------------------------------------------------------------------------
    import paircluster
    import treeutils
    if label == 'data':
        lp_infos = paircluster.read_paired_dir(args.gcreplay_dir)
        with open('%s/selection-metrics.yaml'%args.gcreplay_dir) as sfile:
            smfos = json.load(sfile)
        smdict = {':'.join(u+'-igh' for u in s['unique_ids']) : s for s in smfos}
        antn_pairs = paircluster.get_all_antn_pairs(lp_infos)
        dendro_trees, affy_vals = [], {}
        for h_atn, l_atn in antn_pairs:
            if len(h_atn['unique_ids']) < args.min_seqs_per_gc:
                continue
            gcstr = h_atn['unique_ids'][0].split('-')[0]
            assert gcstr[:2] == 'gc'
            gcn = int(gcstr[2:])
            if gcn not in args.GCs:
                continue
            sfo = smdict.get(':'.join(h_atn['unique_ids']))
            dtree = treeutils.get_dendro_tree(treestr=sfo['lb']['tree'])
            dendro_trees.append(dtree)
            for node in dtree.preorder_node_iter():
                affy_vals[node.taxon.label] = utils.per_seq_val(h_atn, 'affinities', node.taxon.label+'-igh')
        plotvals = get_plotvals(label, dendro_trees, affy_vals)
    elif label == 'simu':
        with open('%s/meta.json'%args.simu_dir) as mfile:
            mfos = json.load(mfile)
        dendro_trees = [treeutils.get_dendro_tree(treestr=s) for s in treeutils.get_treestrs_from_file('%s/trees.nwk'%args.simu_dir)]
        plotvals = get_plotvals(label, dendro_trees, {u : mfos[u]['affinity'] for u in mfos})
    else:
        assert False
    hists = {}
    for pkey, pvals in plotvals.items():
        htmp = Hist(xmin=-10, xmax=5, n_bins=30, value_list=pvals, title=label, xtitle='%s affinity'%pkey)
        htmp.title += ' (%d nodes in %d trees)' % (len(pvals), len(dendro_trees))
        hists['distr-%s-affinity'%pkey] = htmp
    return hists

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
    fn = plotting.draw_no_root(None, plotdir=plotdir, plotname='%s-%s'%(hname, abtype), more_hists=hists, log='y' if abtype=='abundances' else '', xtitle=hists[0].xtitle, ytitle=ytitle,
                               bounds=xbounds, ybounds=ybounds, xticks=xticks, yticks=yticks, yticklabels=yticklabels, errors=hname!='max', square_bins=hname=='max', linewidths=[4, 3],
                               plottitle='mean distr. over GCs' if 'N seqs in bin' in ytitle else '',  # this is a shitty way to identify the mean_hdistr hists, but best i can come up with atm
                               alphas=[0.6, 0.6], colors=[colors[l] for l in labels], translegend=[-0.65, 0] if 'affinity' in abtype else [-0.2, 0], write_csv=True, hfile_labels=labels)
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
parser.add_argument('--gcreplay-dir', default='/fh/fast/matsen_e/processed-data/partis/taraki-gctree-2021-10/v13/delta_bind_CGG_FVS_additive/all')  # default='%s/projects/gcreplay'%partis_dir)
parser.add_argument('--simu-dir')
parser.add_argument('--outdir')
parser.add_argument('--min-seqs-per-gc', type=int, default=70)
parser.add_argument('--mice', default=[1, 2, 3, 4, 5, 6], help='restrict to these mouse numbers')
parser.add_argument('--GCs', default=[0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 55, 56, 57, 58, 59, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83], help='restrict to these GC numbers (these correspond to --mice, used command in comment: ') # csv -cHK_key_gc:HK_key_mouse /fh/fast/matsen_e/data/taraki-gctree-2021-10/gcreplay/nextflow/results/latest/merged-results/observed-seqs.csv |sort|uniq|grep ' [123456]$'|ap 1|sort|uniq >/tmp/out')
parser.add_argument('--is-simu', action='store_true')
parser.add_argument('--gcdyn-dir', default='%s/work/partis/projects/gcdyn'%os.getenv('HOME'))
parser.add_argument('--max-gc-plots', type=int, default=0, help='only plot individual (per-GC) plots for this  many GCs')
args = parser.parse_args()

dlabels = []
if args.data_dir is not None:
    dlabels.append([args.data_dir, 'data'])
if args.simu_dir is not None:
    dlabels.append([args.simu_dir, 'simu'])

abtypes = ['leaf-affinity', 'internal-affinity', 'abundances', 'hdists', 'max-abdn-shm']
hclists = {t : {'distr' : [], 'max' : []} for t in abtypes}
fnames = [[], [], []]
for idir, tlab in dlabels:
    parse_fastas(idir, tlab, is_simu=tlab=='simu')  # runs abundance.py to parse input fastas into a summary csv format
    utils.prep_dir('%s/plots/%s'%(args.outdir, tlab), wildlings=['*.csv', '*.svg'])
    if any('affinity' in t for t in abtypes):
        lhists = get_affinity_plots(tlab)
        for tk in lhists:
            hn, ntype, astr = tk.split('-')
            hclists[ntype+'-'+astr][hn].append(lhists[tk])
    for abtype in [t for t in abtypes if 'affinity' not in t]:
        lhists = plot_abdn_stuff('%s/plots/%s'%(args.outdir, tlab), tlab, abtype)
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

if len(fnames[0]) > 4:
    fnames = [fnames[0][:4], fnames[0][4:]] + fnames[1:]
plotting.make_html(args.outdir+'/plots', fnames=fnames)
dfn = '%s/diff-vals.yaml' % args.outdir
with open(dfn, 'w') as dfile:
    json.dump(diff_vals, dfile)
