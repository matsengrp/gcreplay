#!/usr/bin/env python

from utils import *
import argparse
import ete3
from historydag import beast_loader
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt
import os

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
    naive_sequence=naive_hk_bcr_nt 
    trees = list(beast_loader.load_beast_trees(args.xml_file, args.nexus_file, reference_sequence=naive_sequence)[0])[1:]
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

    def phenotypes(sequence):
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

    def dendropy2ete(tree, naive="naive@0"):
        assert tree.is_rooted

        # determine youngest tip time
        last_time = 0.0
        for leaf in tree.leaf_node_iter():
            if leaf.taxon:
                time = float(leaf.taxon.label.split("@")[1])
                if time > last_time:
                    last_time = time

        # set node ages
        for node in tree.postorder_node_iter():
            if node.taxon:
                node.age = last_time - float(node.taxon.label.split("@")[1])
            else:
                ages = [child.age + child.edge_length for child in node.child_nodes()]
                assert np.allclose(ages, ages[0]), ages
                node.age = ages[0]

        ete_tree = ete3.Tree(dist=0.0)
        for dendropy_node, ete_node in zip(tree.preorder_node_iter(), ete_tree.traverse(strategy="preorder")):
            ete_node.sequence = dendropy_node.comments[0].strip('&states="')
            ete_node.phenotypes = phenotypes(ete_node.sequence)
            ete_node.age = dendropy_node.age
            ete_node.mutations = []  # for tracking mutations on parent branch
            style = ete3.NodeStyle()
            style["hz_line_width"] = 3
            style["vt_line_width"] = 3
            if dendropy_node.is_leaf():
                # ete_node.name = dendropy_node.taxon.label
                if dendropy_node.taxon.label == naive:
                    style["fgcolor"] = "red"
                    naive_node = ete_node
                    style["size"] = 10
                else:
                    style["fgcolor"] = "black"
                    style["size"] = 0
            else:
                style["size"] = 0
            for child in dendropy_node.child_nodes():
                ete_child = ete3.Tree(name=child.label, dist=child.edge_length)
                ete_node.add_child(ete_child)
            ete_node.set_style(style)
            if len(dendropy_node.comments) == 2:
                mutation_history = eval(dendropy_node.comments[1].strip('&history_all=').replace("{", "[").replace("}", "]").replace("A", "'A'").replace("C", "'C'").replace("G", "'G'").replace("T", "'T'"))
                mutation_history.sort(key=lambda x: x[1], reverse=True)  # sort by time
                hamming_distance = sum(x != y for x, y in zip(ete_node.sequence, ete_node.up.sequence))
                assert hamming_distance <= len(set(pos for pos, time, anc, der in mutation_history)), (hamming_distance, mutation_history)
                parent_sequence = list(ete_node.up.sequence)  # need this to track sequence changes in temporal order
                for pos, age, ancestral, derived in mutation_history:
                    pos -= 1  # BEAST is 1-based
                    assert dendropy_node.age < age < dendropy_node.parent_node.age
                    assert parent_sequence[pos] == ancestral, (pos, parent_sequence[pos], ancestral, derived)
                    parent_sequence[pos] = derived  # update parent sequence with this mutation event
                    ete_node.mutations.append((pos, age, ancestral, derived, phenotypes("".join(parent_sequence))))

        if not naive_node.up.is_root():
            print("WARNING: naive sequence is not an outgroup")

        return ete_tree

    def phenotype_slice(tree, age):
        def is_leaf_fn(node):
            return node.age <= age and (node.up is None or node.up.age > age)

        slice_population = []
        for node in tree.iter_leaves(is_leaf_fn=is_leaf_fn):
            if node.age == age:
                slice_population.append(node.phenotypes)
            else:
                current_phenotypes = node.up.phenotypes
                for pos, mutation_age, ancestral, derived, phenotypes in node.mutations:
                    if mutation_age < age:
                        break
                    current_phenotypes = phenotypes
                slice_population.append(current_phenotypes)
        result = pd.concat(slice_population, axis=1).T
        return result

    tree_start_idx=int(len(trees)*args.burn_frac)
    slices = []
    nodes = []
    for tree_idx, tree in enumerate(trees[tree_start_idx:], tree_start_idx):
        ete_tree = dendropy2ete(tree)
        pickle.dump(ete_tree, open(f"{args.outdir}/tree-{tree_idx}.pkl", "wb"))

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
                for pos, age, ancestral, derived, mutation_phenotypes in node.mutations:
                    nodes.append([tree_idx, age, mutation_phenotypes.delta_bind_CGG, mutation_phenotypes.delta_expr])
            nodes.append([tree_idx, node.age, node.phenotypes.delta_bind_CGG, node.phenotypes.delta_expr])


    slice_df = pd.concat(slices, axis=1).T
    slice_df["time"] = slice_df.age.max() - slice_df.age
    slice_df.to_csv(f"{args.outdir}/slice_df.csv")

    node_df = pd.DataFrame(nodes, columns=["tree_idx", "age", "delta_bind_CGG", "delta_expr"])
    node_df["time (days post-immunization)"] = pd.cut(
            node_df.age.max() - node_df.age, bins=100
    ).apply(lambda x: 0 if x.left < 0 else x.right)
    node_df.to_csv(f"{args.outdir}/node_df.csv")
