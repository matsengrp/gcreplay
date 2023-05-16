#!/usr/bin/env python

from utils import *
import beast_utils
import argparse
import ete3
from historydag import beast_loader
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt
import os
from functools import partial

parser = argparse.ArgumentParser(description='A script that takes in string and file arguments')

# Adding string arguments
parser.add_argument('--xml_file', type=str, help='')
parser.add_argument('--nexus_file', type=str, help='')
parser.add_argument('--dms_df', type=str, help='')
parser.add_argument('--pos_df', type=str, help='')
parser.add_argument('--outdir', type=str, help='')
parser.add_argument('--burn_frac', type=float, required=False, default=0.9, help='')
parser.add_argument('--igk_idx', type=int, required=False, default=336, help='')
parser.add_argument('--igh_frame', type=int, required=False, default=1, help='')
parser.add_argument('--igk_frame', type=int, required=False, default=1, help='')


if __name__ == '__main__':

    args = parser.parse_args()

    if not os.path.exists(args.outdir): os.mkdir(args.outdir)
    if not os.path.exists(f"{args.outdir}/trees"): os.mkdir(f"{args.outdir}/trees")
    naive_sequence=naive_hk_bcr_nt 


    # Gather and wrangle DMS data for phenotyping function
    dms_df = pd.read_csv(args.dms_df, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype()))
    # remove linker sites
    dms_df = dms_df[dms_df.chain != "link"]
    dms_df["WT"] = dms_df.wildtype == dms_df.mutant
    assert dms_df.position_IMGT.max() < 1000
    dms_df["site"] = [
        f"{chain}-{str(pos).zfill(3)}" for chain, pos in zip(dms_df.chain, dms_df.position_IMGT)
    ]

    pos_df = pd.read_csv(args.pos_df, dtype=dict(site=pd.Int16Dtype()), index_col="site_scFv")

    # position maps for heavy and light chain that can be used in gctree render
    igh_pos_map = pos_df.loc[pos_df.chain
                             == "H", "site"].reset_index(drop=True)
    igk_pos_map = pos_df.loc[pos_df.chain
                             == "L", "site"].reset_index(drop=True)

    igk_idx, igh_frame, igk_frame = args.igk_idx, args.igh_frame, args.igk_frame

    naive_igh_aa = aa(naive_hk_bcr_nt[:igk_idx], igh_frame)
    naive_igk_aa = aa(naive_hk_bcr_nt[igk_idx:], igk_frame)

    def phenotypes(
        igk_idx,
        naive_igh_aa,
        naive_igk_aa,
        igh_pos_map,
        igk_pos_map,
        sequence
    ):
        igh_aa = aa(sequence[:igk_idx], igh_frame)
        igk_aa = aa(sequence[igk_idx:], igk_frame)
        igh_mutations = mutations(naive_igh_aa, igh_aa, igh_pos_map, "(H)")
        igk_mutations = mutations(naive_igk_aa, igk_aa, igk_pos_map, "(L)")
        all_mutations = igh_mutations + igk_mutations
        fields = ["delta_bind_CGG", "delta_expr"]
        has_stop = any("*" in mutation for mutation in all_mutations)
        if has_stop:
            return pd.Series([dms_df[field].min() for field in fields], index=fields)
        return dms_df.loc[all_mutations, fields].sum(0)


    def phenotype_slice(tree, age):
        def is_leaf_fn(node):
            return node.age <= age and (node.up is None or node.up.age > age)

        slice_population = []
        for node in tree.iter_leaves(is_leaf_fn=is_leaf_fn):
            if node.age == age:
                slice_population.append(node.phenotypes)
            else:
                current_phenotypes = node.up.phenotypes
                for mut_event in node.mutations:

                    if mut_event.age < age:
                        break
                    current_phenotypes = mut_event.phenotypes

                slice_population.append(current_phenotypes)

        result = pd.concat(slice_population, axis=1).T
        return result

    tree_start_idx, ete_trees = beast_utils.beast_trees_as_ete(
        args.xml_file, 
        args.nexus_file, 
        naive_sequence,
        partial(phenotypes, igk_idx, naive_igh_aa, naive_igk_aa, igh_pos_map, igk_pos_map),
        burn_frac=args.burn_frac
    )

    slices = []
    nodes = []
    for tree_idx, ete_tree in enumerate(ete_trees, tree_start_idx):
        pickle.dump(ete_tree, open(f"{args.outdir}/trees/{tree_idx}.pkl", "wb"))

        # phenotype slices through tree time
        for age in np.linspace(0, 20, 20):
            slice = phenotype_slice(ete_tree, age)
            slice_summary = slice.mean()
            slice_summary["age"] = age
            slice_summary["lineages"] = len(slice)
            slice_summary["tree_idx"] = tree_idx
            slices.append(slice_summary)

        # aggregated node data
        for node in ete_tree.traverse():
            if node.mutations:
                for mut_event in node.mutations:
                    nodes.append(
                        [
                            tree_idx, 
                            mut_event.age, 
                            mut_event.phenotypes.delta_bind_CGG, 
                            mut_event.phenotypes.delta_expr
                        ]
                    )

            nodes.append([tree_idx, node.age, node.phenotypes.delta_bind_CGG, node.phenotypes.delta_expr])

    slice_df = pd.concat(slices, axis=1).T
    slice_df["time"] = slice_df.age.max() - slice_df.age
    slice_df["gc"] = "-".join(args.outdir.split("-")[1:])
    slice_df.to_csv(f"{args.outdir}/slice_df.csv", index=False)

    node_df = pd.DataFrame(nodes, columns=["tree_idx", "age", "delta_bind_CGG", "delta_expr"])
    node_df["time (days post-immunization)"] = pd.cut(
            node_df.age.max() - node_df.age, bins=100
    ).apply(lambda x: 0 if x.left < 0 else x.right)
    node_df.to_csv(f"{args.outdir}/node_df.csv")
