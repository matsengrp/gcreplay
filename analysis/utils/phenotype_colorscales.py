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
    cmap.set_bad(color="#ffd400")


class affinity_trees:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")


class affinity_trees_grey_centered:
    colors = ["#A94E35", "#F48365", "#D3D3D3", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")


class expression_dms:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-2.25, vmax=0.25)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")


class expression_trees:
    colors = ["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=0.5)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")

class expression_trees_grey_centered:
    colors = ["#A94E35", "#F48365", "#D3D3D3", "#7378B9", "#383C6C"]
    norm = mc.Normalize(vmin=-3, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")


class psr_dms:
    colors = list(reversed(["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]))
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=2.25)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")


class psr_trees:
    colors = list(reversed(["#A94E35", "#F48365", "#FFFFFF", "#7378B9", "#383C6C"]))
    norm = mc.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=3)
    positions = np.linspace(0, 1, len(colors))
    cmap = mc.LinearSegmentedColormap.from_list(
        "custom_cmap", list(zip(positions, colors))
    )
    cmap.set_bad(color="#ffd400")
