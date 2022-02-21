import pandas as pd
import pickle
from Bio.Seq import Seq

# ARGUMENTS
# ---------
# Path to DMS repo
dms_path = "../../Ab-CGGnaive_DMS"
# path to a gctree inference pickle file
gctree_file = (
    "../MIGRATED_FROM_Ab-CGGnaive_DMS/data/PR1.2c/gc38HK-2/gctree.out.inference.1.p"
)
# output directory
output_dir = "."
# chain reading frames
igh_frame = 1
igk_frame = 1
# IgH+IgK concatenation position
igk_idx = 336
# LBI parameters
tau = 1
tau0 = 0.1
# ---------

# DMS single mutant scores
dms_df = pd.read_csv(
    f"{dms_path}/results/final_variant_scores/final_variant_scores.csv",
    index_col="mutation",
    dtype=dict(position_IMGT=pd.Int16Dtype()),
)
# remove linker sites
dms_df = dms_df[dms_df.chain != "link"]

# position maps for scFv
pos_df = pd.read_csv(
    f"{dms_path}/data/CGGnaive_sites.csv",
    dtype=dict(site=pd.Int16Dtype()),
    index_col="site_scFv",
)
# position maps for heavy and light chain that can be used in gctree render
igh_pos_map = pos_df.loc[pos_df.chain == "H", "site"].reset_index(drop=True)
igk_pos_map = pos_df.loc[pos_df.chain == "L", "site"].reset_index(drop=True)

# load tree object
with open(gctree_file, "rb") as f:
    tree = pickle.load(f)

# generate LBI and LBR node features
tree.local_branching(tau=tau, tau0=tau0)

# Some convenient functions
def aa(seq, frame):
    """Amino acid translation of nucleotide sequence in frame 1, 2, or 3."""
    return Seq(
        seq[(frame - 1) : (frame - 1 + (3 * ((len(seq) - (frame - 1)) // 3)))]
    ).translate()


def mutations(naive_aa, aa, pos_map, chain_annotation):
    """Amino acid substitutions between two sequences, in DMS coordinates."""
    return [
        f"{aa1}{pos_map[pos]}{chain_annotation}{aa2}"
        for pos, (aa1, aa2) in enumerate(zip(naive_aa, aa))
        if aa1 != aa2
    ]


phenotypes = ["delta_bind", "delta_expr", "delta_psr"]

# generate additive phenotype node features, and rows of node dataframe
dat = []
naive_igh_aa = aa(tree.tree.sequence[:igk_idx], igh_frame)
naive_igk_aa = aa(tree.tree.sequence[igk_idx:], igk_frame)
for node in tree.tree.traverse():
    igh_aa = aa(node.sequence[:igk_idx], igh_frame)
    igk_aa = aa(node.sequence[igk_idx:], igk_frame)
    igh_mutations = mutations(naive_igh_aa, igh_aa, igh_pos_map, "(H)")
    igk_mutations = mutations(naive_igk_aa, igk_aa, igk_pos_map, "(L)")
    all_mutations = igh_mutations + igk_mutations
    node.add_feature("mutations", all_mutations)
    has_stop = any("*" in x for x in all_mutations)
    row = [
        node.name,
        node.up.name if node.up else None,
        node.abundance,
        node.sequence[:igk_idx],
        str(igh_aa),
        node.sequence[igk_idx:],
        str(igk_aa),
        ",".join(igh_mutations),
        ",".join(igk_mutations),
        node.LBI,
        node.LBR,
    ]
    for phenotype in phenotypes:
        node.add_feature(
            phenotype, None if has_stop else dms_df.loc[all_mutations, phenotype].sum()
        )
        row.append(getattr(node, phenotype))
    dat.append(row)

# write dataframe
columns = [
    "name",
    "parent_name",
    "abundance",
    "IgH_nt_sequence",
    "IgH_aa_sequence",
    "IgK_nt_sequence",
    "IgK_aa_sequence",
    "IgH_mutations",
    "IgK_mutations",
    "LBI",
    "LBR",
] + phenotypes
df = pd.DataFrame(dat, columns=columns).set_index("name")
df.to_csv(f"{output_dir}/node_data.csv")

# render tree with colormapped features
for phenotype in phenotypes + ["LBI", "LBR"]:
    if phenotype.startswith("LB"):
        cmap = "plasma"
        vmin = 0
        vmax = 10
    if phenotype.startswith("delta_"):
        cmap = "coolwarm"
        # psr has opposite colormap because increase means it is worse
        if phenotype != "delta_psr":
            cmap += "_r"
        vmin = -2
        vmax = 2
    colormap = tree.feature_colormap(phenotype, cmap=cmap, vmin=vmin, vmax=vmax)
    # need this loop to render to svg and notebook
    rendering = tree.render(
        f"{output_dir}/{phenotype}.svg",
        # scale=None, branch_margin=-7,
        scale=20,
        frame=igh_frame,
        frame2=igk_frame,
        chain_split=igk_idx,
        colormap=colormap,
        position_map=igh_pos_map,
        position_map2=igk_pos_map,
    )

# write the new featurized tree to a pickle file
with open(f"{output_dir}/gctree.p", "wb") as f:
    pickle.dump(tree, f)
