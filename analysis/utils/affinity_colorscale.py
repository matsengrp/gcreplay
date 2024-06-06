"""Color scales for affinity values."""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np
# import pandas as pd
# import seaborn as sns


class dms:
    """Based on Tyler's R code for DMS heatmaps"""
    colors = ["#A94E35", "#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = norm = mc.Normalize(vmin=-3.5, vmax=1)
    positions = np.array([0, 1, 2.25, 3.5, 4, 4.5]) / 4.5
    cmap = mc.LinearSegmentedColormap.from_list("custom_cmap",
                                                list(zip(positions, colors)))
    cmap.set_bad(color="#b3b3b3")


class trees:
    norm = mc.CenteredNorm(vcenter=0, halfrange=3)
    cmap = plt.cm.get_cmap("seismic_r")
    cmap.set_bad(color="#b3b3b3")


# tree_delta_bind = pd.Series([node.delta_bind for tree in trees.values() for node in tree.tree.traverse() if not np.isnan(node.delta_bind)])
# delta_bind_CGG = pd.read_csv("https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/improved-Kd-fitting/tite-seq-modeling/output/final_variant_scores.csv", index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())).delta_bind_CGG
# lb = -3
# ub = 3
# sns.histplot(x=tree_delta_bind, label="trees")
# sns.histplot(x=delta_bind_CGG, label="DMS")
# plt.axvline(lb, color='r')
# plt.axvline(ub, color='r')
# plt.legend()
# plt.show()