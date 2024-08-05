"""Color scales for affinity values."""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np


class affinity_dms:
    """Based on Tyler's R code for DMS heatmaps"""

    colors = ["#A94E35", "#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3.5, vmax=1)
    positions = np.array([0, 1, 2.25, 3.5, 4, 4.5]) / 4.5
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class affinity_trees:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class affinity_trees_grey_middle:
    colors = ["#A94E35", "#F48365", "#808080", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class expression_dms:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-2.25, vmax=0.25)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class expression_trees:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=0.5)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class psr_dms:
    colors = list(reversed(["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]))
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=2.25)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


class psr_trees:
    colors = list(reversed(["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]))
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#b3b3b3")


def main():
    import pandas as pd
    import seaborn as sns
    import pathlib
    import glob
    import pickle

    results = pathlib.Path(__file__).parent / "../../nextflow/results"
    files = glob.glob(f"{results}/gctrees/PR*/gctree.p")
    trees = [pickle.load(open(file, "rb")) for file in files]

    dms = "https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/improved-Kd-fitting/tite-seq-modeling/output/final_variant_scores.csv"

    _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

    sns.histplot(
        pd.read_csv(
            dms, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
        ).delta_bind_CGG,
        label="DMS",
        ax=ax1,
    )
    sns.histplot(
        pd.Series(
            [
                node.delta_bind
                for tree in trees
                for node in tree.tree.traverse()
                if not np.isnan(node.delta_bind)
            ]
        ),
        label="trees",
        zorder=0,
        ax=ax1,
    )
    ax1.legend()

    sns.histplot(
        pd.read_csv(
            dms, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
        ).delta_expr,
        label="DMS",
        ax=ax2,
    )
    sns.histplot(
        pd.Series(
            [
                node.delta_expr
                for tree in trees
                for node in tree.tree.traverse()
                if not np.isnan(node.delta_expr)
            ]
        ),
        label="trees",
        zorder=0,
        ax=ax2,
    )
    ax2.legend()

    sns.histplot(
        pd.read_csv(
            dms, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
        ).delta_psr,
        label="DMS",
        ax=ax3,
    )
    sns.histplot(
        pd.Series(
            [
                node.delta_psr
                for tree in trees
                for node in tree.tree.traverse()
                if not np.isnan(node.delta_psr)
            ]
        ),
        label="trees",
        zorder=0,
        ax=ax3,
    )
    ax3.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
