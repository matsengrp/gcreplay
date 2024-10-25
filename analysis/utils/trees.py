r"""Utilities for loading and working with the tree data."""

import glob
import pickle
import ete3
import gctree
import pandas as pd
import numpy as np
import gctree


def load_trees(
    metadata: pd.DataFrame, results: str, ranking_subdir: str = "default"
) -> dict[tuple[str, str, str], ete3.Tree]:
    r"""Load the trees from the given results directory.

    Args:
        metadata: Metadata DataFrame.
        results: Path to the results directory.
        ranking_subdir: GCtree ranking subdirectory.

    Returns:
        A dictionary mapping (mouse, GC) to the corresponding tree.
    """
    trees = {}
    for idx in metadata.index:
        mouse = metadata.mouse[idx]
        gc = metadata.gc[idx]
        file = glob.glob(
            f"{results}/gctrees/PR*-{mouse}-*-{gc}-*/{ranking_subdir}/gctree.p"
        )
        if not file:
            print("missing tree for", idx)
            continue
        if len(file) > 1:
            raise ValueError(f"multiple trees found for {idx}: {file}")
        file = file[0]
        trees[idx] = pickle.load(open(file, "rb"))
    return trees


def burst_stat(tree: gctree.CollapsedTree, stat: str, tau: float) -> None:
    r"""Calculate the specified node statistic and store it in the tree nodes.

    Args:
        tree: The tree.
        stat: The statistic to calculate.
        tau: The tau parameter for the LB or REI statistic.
    """
    if stat == "REI":
        total_abundance = sum(node.abundance for node in tree.tree.traverse())
        for node in tree.tree.traverse():
            REI = (
                sum(
                    node2.abundance * tau ** node2.get_distance(node)
                    for node2 in node.traverse()
                )
                / total_abundance
            )
            node.REI = np.clip(REI, a_min=0.0, a_max=None)
    elif stat.startswith("LB"):
        tree.local_branching(tau=tau)
    else:
        raise ValueError(f"unknown stat {stat}")
