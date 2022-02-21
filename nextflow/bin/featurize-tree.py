import click
import pandas as pd
import pickle
from Bio.Seq import Seq


def aa(sequence, frame):
    """Amino acid translation of nucleotide sequence in frame 1, 2, or 3."""
    return Seq(
        sequence[(frame - 1) : (frame - 1 + (3 * ((len(sequence) - (frame - 1)) // 3)))]
    ).translate()


def mutations(naive_aa, aa, pos_map, chain_annotation):
    """Amino acid substitutions between two sequences, in IMGT coordinates."""
    return [
        f"{aa1}{pos_map[pos]}{chain_annotation}{aa2}"
        for pos, (aa1, aa2) in enumerate(zip(naive_aa, aa))
        if aa1 != aa2
    ]


@click.command()
@click.argument("gctree_file", type=click.Path(exists=True))
@click.argument("variant_scores", type=click.Path(exists=True))
@click.argument("naive_sites", type=click.Path(exists=True))
@click.option(
    "--igk_idx",
    "-k",
    type=int,
    required=True,
    help="Start index of light chain in concatenated sequence.",
)
@click.option("--igh_frame", type=int, default=1, help="Heavy chain reading frame.")
@click.option("--igk_frame", type=int, default=1, help="Light chain reading frame.")
@click.option("--tau", type=float, default=1.0, help="Bandwidth for LBI kernel.")
@click.option(
    "--tau0",
    type=float,
    default=0.1,
    help="Branch length assigned to zero-mutation branches for LBI computation.",
)
@click.option(
    "--output_dir",
    "-o",
    type=click.Path(exists=True),
    default=".",
    help="Path to write output files.",
)
def main(
    gctree_file,
    variant_scores,
    naive_sites,
    igk_idx,
    igh_frame,
    igk_frame,
    tau,
    tau0,
    output_dir,
):
    """Featurizes a gctree.CollapsedTree object using DMS data, outputing a new
    pickled tree object, a csv of node data, and rendered trees according to features.

    \b
    GCTREE_FILE: Path to pickled gctree.CollapsedTree object
    VARIANT_SCORES: Path to variant scores csv in DMS repository
    NAIVE_SITES: Path to sites csv in DMS repository
    """
    # DMS single mutant scores
    dms_df = pd.read_csv(
        variant_scores, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
    )
    # remove linker sites
    dms_df = dms_df[dms_df.chain != "link"]

    # position maps for scFv
    pos_df = pd.read_csv(
        naive_sites,
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
                phenotype,
                None if has_stop else dms_df.loc[all_mutations, phenotype].sum(),
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
        elif phenotype.startswith("delta_"):
            cmap = "coolwarm"
            # psr has opposite colormap because increase means it is worse
            if phenotype != "delta_psr":
                cmap += "_r"
            vmin = -2
            vmax = 2
        else:
            raise RuntimeError("unknown phenotype")
        colormap = tree.feature_colormap(phenotype, cmap=cmap, vmin=vmin, vmax=vmax)
        # need this loop to render to svg and notebook
        tree.render(
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


if __name__ == "__main__":
    main()
