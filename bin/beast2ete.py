#!/usr/bin/env python

from utils import aa, mutations, naive_hk_bcr_nt
from beast_utils import beast_trees_as_ete

import pickle
import argparse
import pandas as pd
import numpy as np
import os
import re
from functools import partial

# TODO remove
# import ete3
# from historydag import beast_loader
# from Bio.Seq import Seq
# import seaborn as sns
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="A script that takes in string and file arguments"
)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


# Adding string arguments
parser.add_argument("--xml_file", type=str, help="")
parser.add_argument("--nexus_file", type=str, help="")
parser.add_argument("--dms_df", type=str, help="")
parser.add_argument("--pos_df", type=str, help="")
parser.add_argument("--outdir", type=str, help="")
# parser.add_argument('--pkldir', type=str, required=False, default=None)
parser.add_argument(
    "--save_pkl_trees",
    type=str2bool,
    nargs="?",
    const=True,
    default=False,
    help="Save all MCMC trees as pickles.",
)
parser.add_argument("--burn_frac", type=float, required=False, default=0.9, help="")
parser.add_argument("--igk_idx", type=int, required=False, default=336, help="")
parser.add_argument("--igh_frame", type=int, required=False, default=1, help="")
parser.add_argument("--igk_frame", type=int, required=False, default=1, help="")
# TODO add a parameter to save only a single pkl tree
parser.add_argument(
    "--save_single_pkl",
    type=str2bool,
    nargs="?",
    const=True,
    default=True,
    help="Save the last MCMC tree as a pickle.",
)
parser.add_argument(
    "--n_phenotype_slices",
    type=int,
    required=False,
    default=100,
    help="Number of phenotype slices to take between naive and sampling time. Default is 100.",
)


if __name__ == "__main__":

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    if args.save_pkl_trees or args.save_single_pkl:
        if not os.path.exists(f"{args.outdir}/pkl_ete_trees"):
            os.mkdir(f"{args.outdir}/pkl_ete_trees")

    naive_sequence = naive_hk_bcr_nt

    # Gather and wrangle DMS data for phenotyping function
    dms_df = pd.read_csv(
        args.dms_df, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
    )
    # remove linker sites
    dms_df = dms_df[dms_df.chain != "link"]
    dms_df["WT"] = dms_df.wildtype == dms_df.mutant
    assert dms_df.position_IMGT.max() < 1000
    dms_df["site"] = [
        f"{chain}-{str(pos).zfill(3)}"
        for chain, pos in zip(dms_df.chain, dms_df.position_IMGT)
    ]

    pos_df = pd.read_csv(
        args.pos_df, dtype=dict(site=pd.Int16Dtype()), index_col="site_scFv"
    )

    # position maps for heavy and light chain that can be used in gctree render
    igh_pos_map = pos_df.loc[pos_df.chain == "H", "site"].reset_index(drop=True)
    igk_pos_map = pos_df.loc[pos_df.chain == "L", "site"].reset_index(drop=True)

    igk_idx, igh_frame, igk_frame = args.igk_idx, args.igh_frame, args.igk_frame

    naive_igh_aa = aa(naive_hk_bcr_nt[:igk_idx], igh_frame)
    naive_igk_aa = aa(naive_hk_bcr_nt[igk_idx:], igk_frame)

    def phenotypes(
        igk_idx, naive_igh_aa, naive_igk_aa, igh_pos_map, igk_pos_map, sequence
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

    def phenotype_slice(tree, time):
        def is_leaf_fn(node):
            return node.time >= time and (node.up is None or node.up.time < time)

        if time < 0.0 or time > 1.0:
            raise ValueError(f"time {time} is outside tree time range 0, 1")

        slice_population = []
        for node in tree.iter_leaves(is_leaf_fn=is_leaf_fn):
            if node.time == time:
                slice_population.append(node.phenotypes)
            else:
                current_phenotypes = node.up.phenotypes

                for mut_event in node.mutations:
                    if mut_event["time"] < time:
                        current_phenotypes = mut_event["derived_phenotypes"]
                    else:
                        break

                slice_population.append(current_phenotypes)

        result = pd.concat(slice_population, axis=1).T
        return result

    tree_start_idx, ete_trees = beast_trees_as_ete(
        args.xml_file,
        args.nexus_file,
        naive_sequence,
        partial(
            phenotypes, igk_idx, naive_igh_aa, naive_igk_aa, igh_pos_map, igk_pos_map
        ),
        burn_frac=args.burn_frac,
    )

    slices = []
    nodes = []

    for tree_idx, ete_tree in enumerate(ete_trees, tree_start_idx):
        if args.save_pkl_trees:
            pickle.dump(
                ete_tree, open(f"{args.outdir}/pkl_ete_trees/{tree_idx}.pkl", "wb")
            )
        if args.save_single_pkl and tree_idx == tree_start_idx + len(ete_trees) - 1:
            pickle.dump(
                ete_tree, open(f"{args.outdir}/last_sampled_tree.pkl", "wb")
            )

        for time in np.linspace(0, 1, args.n_phenotype_slices):
            slice = phenotype_slice(ete_tree, time)

            slice["time"] = time.round(3)
            slice["tree_idx"] = tree_idx
            assert np.all(slice.index.values == range(len(slice)))
            slice.index.name = "entry_idx"
            slice.reset_index(drop=False, inplace=True)
            slices.append(slice)

    slice_df = pd.concat(slices).reset_index(drop=True)
    slice_df.to_csv(f"{args.outdir}/slice_df.csv", index=False)
