#!/usr/bin/env python2
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import glob
import argparse
import colored_traceback.always
import json
from Bio.SeqUtils.CheckSum import seguid
import itertools
import pandas as pd
import numpy
import collections
import ast

# if you move this script, you'll need to change this method of getting the imports
partis_dir = '%s/work/partis' % os.getenv('HOME') #os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import plotting
import lbplotting
import hutils
from hist import Hist
import treeutils
import paircluster

# ----------------------------------------------------------------------------------------
colors = {'data' : plotting.default_colors[1],
          'simu' : plotting.default_colors[0]}
pltlabels = {'hdists' : 'root-tip dist', 'max-abdn-shm' : 'SHM in most\nabundant seq'}  # 'median SHM of seqs\nw/max abundance'}

# ----------------------------------------------------------------------------------------
def abfn(tlab, abtype='abundances'):
    return '%s/%s/%s.csv' % (args.outdir, tlab, abtype)

# ----------------------------------------------------------------------------------------
def write_abdn_csv(label, all_seqfos):  # summarize abundance (and other) info into csv files, for later reading (kind of weird to separate it like this, it's because this used to be in another script)
    # ----------------------------------------------------------------------------------------
    def process_family(fam_fos):
        n_seqs = len(fam_fos)
        if args.min_seqs_per_gc is not None and n_seqs < args.min_seqs_per_gc:
            counters['too-small'] += 1
            return
        if args.max_seqs_per_gc is not None and n_seqs > args.max_seqs_per_gc:
            init_sizes.append(n_seqs)
            i_to_keep = numpy.random.choice(list(range(n_seqs)), size=args.max_seqs_per_gc)
            fam_fos = [fam_fos[i] for i in i_to_keep]

        # This dictionary will map sequence checksum to the list of squence ids that have that
        # sequence checksum.
        ids_by_checksum = collections.defaultdict(list)
        hdvals, hd_dict = (
            [],
            {},
        )  # list of hamming distance to naive for each sequence (hd_dict is just for max_abdn below)
        sequence_count = 0
        for sfo in fam_fos:
            sequence_count = sequence_count + 1
            ids_by_checksum[seguid(sfo['seq'])].append(sfo['name'])
            hdvals.append(sfo['n_muts'])
            hd_dict[sfo['name']] = hdvals[-1]

        abundance_distribution = collections.defaultdict(int)
        for id_list in ids_by_checksum.values():
            id_count = len(id_list)
            abundance_distribution[id_count] = abundance_distribution[id_count] + 1
        assert sequence_count == sum(k * v for k, v in abundance_distribution.items())
        for amax, max_abdn_idlists in itertools.groupby(
            sorted(ids_by_checksum.values(), key=len, reverse=True), key=len
        ):
            # print(amax, [len(l) for l in max_abdn_idlists])
            break  # just want the first one (max abundance)

        base = hash(''.join(s['seq'] for s in fam_fos))
        abundances[base] = pd.Series(
            abundance_distribution.values(), index=abundance_distribution.keys()
        )

        fdicts["hdists"][base] = hdvals
        fdicts["max-abdn-shm"][base] = [
            int(numpy.median([hd_dict[u] for x in max_abdn_idlists for u in x]))
        ]

    # ----------------------------------------------------------------------------------------
    counters = {'too-small' : 0, 'removed' : 0}
    init_sizes = []
    abundances = {}
    fdicts = {"hdists": {}, "max-abdn-shm": {}}
    for gcn, fam_fos in all_seqfos.items():
        process_family(fam_fos)

    if counters['too-small'] > 0:
        print("    skipped %d files with fewer than %d seqs" % (counters['too-small'], args.min_seqs_per_gc))
    if len(init_sizes) > 0:
        print(
            "    downsampled %d samples to %d from initial sizes: %s"
            % (len(init_sizes), args.max_seqs_per_gc, " ".join(str(s) for s in sorted(init_sizes)))
        )

    print "  writing to %s" % os.path.dirname(abfn(label))
    if not os.path.exists(os.path.dirname(abfn(label))):
        os.makedirs(os.path.dirname(abfn(label)))

    to_write = pd.DataFrame(abundances).fillna(0).astype(int)
    to_write.to_csv(abfn(label, 'abundances'))

    for fstr, fdct in fdicts.items():
        with open(abfn(label, abtype=fstr), 'w') as cfile:
            writer = csv.DictWriter(cfile, ["fbase", "vlist"])
            writer.writeheader()
            for fbase, vlist in fdicts[fstr].items():
                writer.writerow({"fbase": fbase, "vlist": ":".join(str(v) for v in vlist)})

# ----------------------------------------------------------------------------------------
def abdn_hargs(hlist):
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
def plot_abdn_stuff(lhists, plotdir, label, abtype):
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
        hmax.title = '%s (%d trees)' % (label, len(distr_hists))
        hmax.xtitle = 'max abundance in GC'

    # plot mean distribution over GCs
    mean_hdistr = plotting.make_mean_hist(distr_hists)
    mean_hdistr.title = '%s (%d trees)' % (label, len(distr_hists))
    mean_hdistr.xtitle = pltlabels.get(abtype, abtype)
    mean_hdistr.ytitle = 'N seqs in bin\nmean+/-std (over GCs)'

    lhists[abtype] = {'distr' : mean_hdistr, 'max' : hmax}

# ----------------------------------------------------------------------------------------
def read_input_files(label):
    # ----------------------------------------------------------------------------------------
    def nstr(node):
        return 'leaf' if node.is_leaf() else 'internal'
    # ----------------------------------------------------------------------------------------
    def get_simu_affy(label, dendro_trees, affy_vals):
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
    def read_data_file():
        skipped_mice, kept_mice, all_gcs = set(), set(), set()
        print '    reading replay data from %s' % args.gcreplay_dir
        with open('%s/gcreplay/nextflow/results/latest/merged-results/gctree-node-data.csv'%args.gcreplay_dir) as cfile:
            reader = csv.DictReader(cfile)
            for line in reader:
                gcn = int(line['HK_key_gc'])
                all_gcs.add(gcn)
                if int(line['HK_key_mouse']) not in args.mice:
                    skipped_mice.add(int(line['HK_key_mouse']))
                    continue
                kept_mice.add(int(line['HK_key_mouse']))
                if gcn not in all_seqfos:
                    all_seqfos[gcn] = []
                affinity = line['delta_bind_CGG_FVS_additive']
                hdist = utils.hamming_distance(args.naive_seq, line['IgH_nt_sequence']+line['IgK_nt_sequence'])  # this is different to nmuts in the next line, not yet sure why (UPDATE: one is nuc/other is AA, docs/column names maybe will be changed in gcreplay)
                # nmuts = int(line['n_mutations_HC']) + int(line['n_mutations_LC'])
                # print '  %2d  %2d  %s' % (hdist, nmuts, utils.color('red', '<--') if hdist!=nmuts else '')
                # print '      %s' % utils.color_mutants(NAIVE_SEQUENCE, line['IgH_nt_sequence']+line['IgK_nt_sequence'])
                for iseq in range(int(line['abundance'])):
                    all_seqfos[gcn].append({'name' : 'gc%d-%s-%d' % (int(line['HK_key_gc']), line['name'], iseq),
                                            'seq' : line['IgH_nt_sequence']+line['IgK_nt_sequence'],
                                            # 'n_muts' : nmuts,
                                            'n_muts' : hdist,
                                            'affinity' : None if affinity == '' else float(affinity),
                                            })
        print '    kept %d / %d GCs from  %d / %d mice: %s' % (len(all_seqfos), len(all_gcs), len(kept_mice), len(kept_mice) + len(skipped_mice), ' '.join(str(m) for m in kept_mice))
    # ----------------------------------------------------------------------------------------
    all_seqfos = collections.OrderedDict()
    plotvals = {k : [] for k in ['leaf', 'internal']}
    n_missing, n_tot = {'internal' : [], 'leaf' : []}, {'internal' : [], 'leaf' : []}  # per-seq (not per-gc) counts
    if label == 'data':
        read_data_file()  # read seqs plus affinity and mutation info from csv file (still have to read trees below to get leaf/internal info)
        n_too_small = 0
        for gcn in all_seqfos:
            if len(all_seqfos[gcn]) < args.min_seqs_per_gc:  # NOTE not subsampling with args.max_seqs_per_gc as in write_abdn_csv() (can't be bothered, it seems more important for abundance stuff anyway)
                n_too_small += 1
                continue
            gctree_dir = utils.get_single_entry(glob.glob('%s/gcreplay/nextflow/results/latest/gctrees/PR*-%d-GC'%(args.gcreplay_dir, gcn)))
            dtree = treeutils.get_dendro_tree(treefname='%s/gctree.inference.1.nk'%gctree_dir)
            nodefo = {n.taxon.label : n for n in dtree.preorder_node_iter()}
            for sfo in all_seqfos[gcn]:
                gcstr, nname, iseq = sfo['name'].split('-')
                n_tot[nstr(nodefo[nname])].append(nname)
                if sfo['affinity'] is None:
                    n_missing[nstr(nodefo[nname])].append(nname)
                    continue
                plotvals[nstr(nodefo[nname])].append(sfo['affinity'])
        if sum(len(l) for l in n_missing.values()) > 0:
            print '      %s missing/none affinity values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal']))
        if n_too_small > 0:
            print '    skipped %d / %d gcs with fewer than %d seqs' % (n_too_small, len(all_seqfos), args.min_seqs_per_gc)
        n_trees = len(all_seqfos) - n_too_small
    elif label == 'simu':
        with open('%s/meta.json'%args.simu_dir) as mfile:
            mfos = json.load(mfile)
        tmp_seqfos = utils.read_fastx('%s/seqs.fasta'%args.simu_dir)
        for sfo in tmp_seqfos:
            if 'naive' in sfo['name']:
                continue
            sfo['n_muts'] = mfos[sfo['name']]['n_muts']
            gcn, nname = sfo['name'].split('-')
            if gcn not in all_seqfos:
                all_seqfos[gcn] = []
            all_seqfos[gcn].append(sfo)
        dendro_trees = [treeutils.get_dendro_tree(treestr=s) for s in treeutils.get_treestrs_from_file('%s/trees.nwk'%args.simu_dir)]
        plotvals = get_simu_affy(label, dendro_trees, {u : mfos[u]['affinity'] for u in mfos})
        n_trees = len(dendro_trees)
    else:
        assert False

    write_abdn_csv(label, all_seqfos)

    hists = {}
    for pkey, pvals in plotvals.items():
        htmp = Hist(xmin=-15, xmax=7, n_bins=30, value_list=pvals, title=label, xtitle='%s affinity'%pkey)
        htmp.title += ' (%d nodes in %d trees)' % (len(pvals), n_trees)
        hists['%s-affinity'%pkey] = {'distr' : htmp}
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
    ytitle = hists[0].ytitle
    if args.normalize:
        for htmp in hists:
            htmp.normalize()
        if 'fraction of' in hists[0].ytitle:
            ytitle = 'fraction of total'
    xbounds, ybounds, xticks, yticks, yticklabels = abdn_hargs(hists) if abtype=='abundances' else (None, None, None, None, None)  # seems not to need this? hutils.multi_hist_filled_bin_xbounds(hists)
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
parser.add_argument('--gcreplay-dir', default='/fh/fast/matsen_e/data/taraki-gctree-2021-10')
parser.add_argument('--simu-dir')
parser.add_argument('--outdir')
parser.add_argument('--min-seqs-per-gc', type=int, default=70)
parser.add_argument('--max-seqs-per-gc', type=int, default=70)
parser.add_argument('--mice', default=[1, 2, 3, 4, 5, 6], help='restrict to these mouse numbers')
parser.add_argument('--GCs', default=[0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 55, 56, 57, 58, 59, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83], help='restrict to these GC numbers (these correspond to --mice, used command in comment: ') # csv -cHK_key_gc:HK_key_mouse /fh/fast/matsen_e/data/taraki-gctree-2021-10/gcreplay/nextflow/results/latest/merged-results/observed-seqs.csv |sort|uniq|grep ' [123456]$'|ap 1|sort|uniq >/tmp/out')
parser.add_argument('--is-simu', action='store_true')
parser.add_argument('--gcdyn-dir', default='%s/work/partis/projects/gcdyn'%os.getenv('HOME'))
parser.add_argument('--max-gc-plots', type=int, default=0, help='only plot individual (per-GC) plots for this  many GCs')
parser.add_argument('--normalize', action='store_true')
parser.add_argument('--naive-seq', default="GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAGACATTGTGATGACTCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACTAATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCTACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCTCTCACGTTCGGCTCGGGGACTAAGCTAGAAATAAAA")
parser.add_argument(
    "--random-seed", type=int, default=1, help="random seed for subsampling"
)
args = parser.parse_args()

numpy.random.seed(args.random_seed)

dlabels = []
if args.data_dir is not None:
    dlabels.append([args.data_dir, 'data'])
if args.simu_dir is not None:
    dlabels.append([args.simu_dir, 'simu'])

abtypes = ['leaf-affinity', 'internal-affinity', 'abundances', 'hdists', 'max-abdn-shm']
hclists = {t : {'distr' : [], 'max' : []} for t in abtypes}
fnames = [[], [], []]
for idir, tlab in dlabels:
    utils.prep_dir('%s/plots/%s'%(args.outdir, tlab), wildlings=['*.csv', '*.svg'])
    lhists = read_input_files(tlab)  # puts affinity hists in lhists
    for abtype in abtypes:
        if 'affinity' not in abtype:
            plot_abdn_stuff(lhists, '%s/plots/%s'%(args.outdir, tlab), tlab, abtype)  # adds abdn hists to lhists
        for hn in lhists[abtype]:
            hclists[abtype][hn].append(lhists[abtype][hn])

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
