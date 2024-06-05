"""Color scale for affinity values, based on Tyler's R code for DMS heatmaps."""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np

colors = ["#A94E35", "#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
positions = np.array([0, 1, 2.25, 3.5, 4, 4.5]) / 4.5
cmap = mc.LinearSegmentedColormap.from_list("custom_cmap",
                                            list(zip(positions, colors)))
cmap.set_bad(color="#b3b3b3")
norm = mc.Normalize(vmin=-3.5, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])