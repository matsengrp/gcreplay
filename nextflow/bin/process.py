#!/usr/bin/env python
import copy
import sys
import csv
import os
import numpy
import json
import glob
import colored_traceback.always
import re
import argparse

# NOTE run from inside main partis dir
partis_dir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__)).replace('/datascripts/meta/goo-dengue-10x', '')
sys.path.insert(1, partis_dir + '/python')
import utils
import clusterpath
import treeutils
import dendropy

parser = argparse.ArgumentParser()
parser.add_argument('--version', default='v1')  # -single-muts
parser.add_argument('--n-max-aa-muts', type=int, help='if set, either don\'t set kd values (default) or actually skip (if --skip-too-many-mut-seqs is set) abs with more than this many aa mutations (sum over h and l).')
parser.add_argument('--skip-too-many-mut-seqs', action='store_true')
args = parser.parse_args()

kdkeys = ['deltaKd', 'deltaExpression', 'deltaPSR']
outkeys = ['name'] + kdkeys + ['multiplicity']

bdr = '/fh/fast/matsen_e/data/taraki-gctree-2021-10'
rdr = '%s/Ab-CGGnaive_DMS/replay' % bdr
baseoutdir = '%s/processed-data/%s' % (bdr, args.version)
datadirs = ['%s/data/%s'%(rdr, d) for d in ['211015 PR1.1', 'PR1.2c', 'PR1.3', 'PR1.4']]
all_gcdirs = {os.path.basename(gd) : gd for d in datadirs for gd in glob.glob('%s/gc*'%d)}
csvfname = '%s/data/211105PR1_total_select.csv' % rdr
def gctfn(gcn, suffix='.fasta'):
    return '%s/gctree.out.inference.1%s' % (all_gcdirs[gcn], suffix)
def abfn(gcn):
    return '%s/abundances.csv' % all_gcdirs[gcn]

# ----------------------------------------------------------------------------------------
def get_gc_num(gcd):
    gstr = utils.get_single_entry(re.findall('gc[0-9][0-9]*HK', gcd))
    return gstr.lstrip('gc').rstrip('HK')
# ----------------------------------------------------------------------------------------
def get_gc_name(gnum):
    gnames = [g for g in all_gcdirs if 'gc%sHK'%gnum in g]
    if len(gnames) == 0:
        return None
    elif len(gnames) > 1:
        raise Exception('no unique match for %s among %s' % (gnum, ' '.join(all_gcdirs)))
    else:
        return gnames[0]
# ----------------------------------------------------------------------------------------
def split_fstr(fstr):  # it's the two-line str for writing to fasta
    assert '\n' in fstr
    ustr, tmp_seq = fstr.split('\n')
    return tmp_seq
# ----------------------------------------------------------------------------------------
def split_seq(gnum, seq):
    hlen = seq_lens[gnum]['H']
    return seq[:hlen], seq[hlen:]
# ----------------------------------------------------------------------------------------
def addlen(line, gcn, lstr):
    if gcn not in seq_lens:
        seq_lens[gcn] = {}
    if lstr not in seq_lens[gcn]:
        seq_lens[gcn][lstr] = len(split_fstr(line['partis_%s'%lstr]))
# ----------------------------------------------------------------------------------------
def getseq(line, gcn, lstr):
    return split_fstr(line['partis_%s'%lstr])
# ----------------------------------------------------------------------------------------
def nfcn(g, n):
    return 'gc%s-%s' % (g, n)
# ----------------------------------------------------------------------------------------
def add_to_asfos(gnum, sfo):
    asfo = copy.deepcopy(sfo)
    asfo['name'] = nfcn(gnum, sfo['name'])
    all_seqfos.append(asfo)
    if asfo['name'] not in clusters[gnum]:  # we add the h and k seqs with separate calls to this fcn, but they have the same uid at this point
        clusters[gnum].append(asfo['name'])

# ----------------------------------------------------------------------------------------
print '  reading naive seqs + trees from gctree output dirs'
naive_seqs, gtrees, multiplicities = {}, {}, {}
for gcn, gcd in all_gcdirs.items():
    tuple_info = []
    tmpfos = utils.read_fastx(gctfn(gcn), look_for_tuples=True, tuple_info=tuple_info)
    nfos = [s for s in tmpfos if s['name']=='naive']
    if len(nfos) == 0:
        print '  %s no naive seq for %s in %s' % (utils.wrnstr(), gcn, gctfn(gcn))
        continue
    elif len(nfos) > 1:
        print '  %s multiple (%d) naive seqs for %s in %s' % (utils.wrnstr(), len(nfos), gcn, gctfn(gcn))
    nfo = nfos[0]
    gnum = get_gc_num(gcn)
    naive_seqs[gnum] = {'name' : 'naive', 'seq' : nfo['seq'], 'multiplicity' : 1}
    naive_seqs[gnum].update({k : 0 for k in kdkeys})

    # read tree
    dtree = treeutils.get_dendro_tree(treefname=gctfn(gcn, suffix='.nk'))
    dtree.scale_edges(1. / 400)  # NOTE *don't* use the seqs in tmpfos (or naive_seqs), since they're concat'd h+l
    treeutils.translate_labels(dtree, [(dtree.seed_node.taxon.label, 'naive')])
    for nids in tuple_info:  # have to add zero length branches for the little fuckers as well
        tnodes = [dtree.find_node_with_taxon_label(u) for u in nids]
        tnode = utils.get_single_entry([n for n in tnodes if n is not None])
        tns = dtree.taxon_namespace
        for uid in [u for u in nids if u!=tnode.taxon.label]:
            tns.add_taxon(dendropy.Taxon(uid))
            _ = tnode.new_child(taxon=tns.get_taxon(uid), edge_length=1./400)  # ok, maybe not zero. Whatever.
    gtrees[gnum] = dtree

    # read abundances
    if gnum not in multiplicities:
        multiplicities[gnum] = {}
    with open(abfn(gcn)) as afile:
        reader = csv.DictReader(afile, fieldnames=('name', 'abundance'))
        for line in reader:
            multiplicities[gnum][line['name']] = max(1, int(line['abundance']))  # increase 0s (inferred ancestors) to 1
missing_nseqs = set([get_gc_num(d) for d in all_gcdirs.keys()]) - set(naive_seqs)
if len(missing_nseqs) > 0:
    print '    %s missing %d / %d naive seqs (for %s)'  % (utils.wrnstr(), len(missing_nseqs), len(all_gcdirs), ' '.join(missing_nseqs))

# ----------------------------------------------------------------------------------------
# n_skipped, na_strs = 0, set()
n_total, n_multi_mut_skipped, n_dup_skipped = 0, 0, 0
all_seqfos, clusters = [], {}
gc_seqfos, allseqs = {}, {}
seq_lens = {}
bad_names = set()
print '  reading seqs from %s' % csvfname
with open(csvfname) as cfile:
    reader = csv.DictReader(cfile)
    for line in reader:
        gnum = line['GC_num']
        gcn = get_gc_name(gnum)
        if gcn is None:
            bad_names.add(gnum)
            continue
        n_total += 1
        if gnum not in gc_seqfos:
            gc_seqfos[gnum] = []
            allseqs[gnum] = set()
            clusters[gnum] = []
        hk_sfos = []
        for lstr in 'HK':
            addlen(line, gnum, lstr)
            if 'seq' in line['idmap'] and line['idmap'] not in multiplicities[gnum]:
                print '    %s no abundance info for observed seq %s from %s: %s' % (utils.wrnstr(), line['idmap'], gnum, abfn(gcn))
            sfo = {'name' : line['idmap'], 'seq' : getseq(line, gnum, lstr).upper(), 'multiplicity' : multiplicities[gnum].get(line['idmap'], 1)}
            sfo.update({k : line[k] for k in kdkeys})
            hk_sfos.append(sfo)
        def sk(s): return s['name'] + s['seq']
        if any(sk(s) in allseqs[gnum] for s in hk_sfos):  # the number of times we see this is maybe the multiplicity, but sometimes seems different? (if we see either, they should both be in there)
            n_dup_skipped += 1
            continue
        if args.n_max_aa_muts is not None:
            nseq, mseq = naive_seqs[gnum]['seq'], ''.join(s['seq'] for s in hk_sfos)
            # utils.color_mutants(nseq, mseq, print_result=True)
            # utils.color_mutants(utils.ltranslate(nseq), utils.ltranslate(mseq), print_result=True, amino_acid=True)
            hdist = utils.hamming_distance(utils.ltranslate(nseq), utils.ltranslate(mseq), amino_acid=True)
            if hdist > args.n_max_aa_muts:
                n_multi_mut_skipped += 1
                if args.skip_too_many_mut_seqs:
                    continue
                else:
                    for sfo in hk_sfos:
                        for tk in kdkeys:
                            sfo[tk] = None
        for sfo in hk_sfos:
            gc_seqfos[gnum].append(sfo)
            allseqs[gnum].add(sk(sfo))
            add_to_asfos(gnum, sfo)

# if n_skipped > 0:
#     print '  %s skipped %d with unknown fasta strs %s' % (utils.color('yellow', 'warning'), n_skipped, ' '.join(na_strs))
if n_dup_skipped > 0:
    print '    skipped %d / %d lines with duplicate uids+seqs' % (n_dup_skipped, n_total)
if n_multi_mut_skipped > 0:
    tstr = 'skipped' if args.skip_too_many_mut_seqs else 'removed kd info for'
    print '    --n-max-aa-muts: %s %d / %d abs with more than %s aa mutation%s' % (tstr, n_multi_mut_skipped, n_total, args.n_max_aa_muts, utils.plural(args.n_max_aa_muts))
if len(bad_names) > 0:
    print '  %s no gc names correspond to GC_num values %s' % (utils.wrnstr(), ' '.join(bad_names))

# add naive seqs
unused_nfos = naive_seqs.keys()
for gnum in gc_seqfos:
    if gnum in naive_seqs:
        hnfo, lnfo = [copy.deepcopy(naive_seqs[gnum]) for _ in range(2)]
        hnfo['seq'], lnfo['seq'] = split_seq(gnum, naive_seqs[gnum]['seq'])
        gc_seqfos[gnum] = [hnfo, lnfo] + gc_seqfos[gnum]
        unused_nfos.remove(gnum)
        add_to_asfos(gnum, hnfo)
        add_to_asfos(gnum, lnfo)
    else:
        print '  %s no naive seq for %s' % (utils.wrnstr(), gnum)
if len(unused_nfos) > 0:
    print '  didn\'t use %d naive seqs: %s' % (len(unused_nfos), ' '.join(unused_nfos))

# rename nodes in trees for merged file
all_trees = []
for gnum, dtree in gtrees.items():
    ntree = copy.deepcopy(dtree)
    treeutils.translate_labels(ntree, [(n.taxon.label, nfcn(gnum, n.taxon.label)) for n in ntree.preorder_node_iter()])
    all_trees.append(ntree)

print '  writing seqfos to %s' % baseoutdir
utils.write_fasta('%s/all-seqs.fa'%baseoutdir, all_seqfos)
utils.write_only_partition('%s/partition-only.yaml'%baseoutdir, clusters.values())
with open('%s/all-trees.nwk'%baseoutdir, 'w') as tfile:
    for ntree in all_trees:
        tfile.write(treeutils.as_str(ntree) + '\n')
with open('%s/kdvals.csv'%baseoutdir, 'w') as mfile:
    writer = csv.DictWriter(mfile, outkeys)
    writer.writeheader()
    for sfo in all_seqfos:
        writer.writerow({k : sfo[k] for k in outkeys})

for gnum in gc_seqfos:
    ofn = '%s/gc%s/seqs.fa' % (baseoutdir, gnum)
    utils.write_fasta(ofn, gc_seqfos[gnum])
    ofn = '%s/gc%s/kdvals.csv' % (baseoutdir, gnum)
    with open('%s/gc%s/tree.nwk'%(baseoutdir, gnum), 'w') as tfile:
        tfile.write(treeutils.as_str(gtrees[gnum]))
    with open(ofn, 'w') as mfile:
        writer = csv.DictWriter(mfile, outkeys)
        writer.writeheader()
        for sfo in gc_seqfos[gnum]:
            writer.writerow({k : sfo[k] for k in outkeys})
