#!/usr/bin/env python
"""
@file: gctree-tools.py
"""

import pickle
from collections import defaultdict

import numpy as np
import pandas as pd
import click
from click import group

from utils import aa, nt_mutations, mutations, naive_hk_bcr_nt, final_HK_col_order
from trees import burst_stat


#################################
# CLI
#################################


# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Welcome to the gcreplay-tools CLI!

    Here we present a few useful utilities for the
    processing and analysis of BCR's extracted from
    a gcreplay experiment. for any one the sub-commands,
    listed below, you may use gcreplay-tools <command>
    --help for more information on the parameters of interest.
    """
    pass

@cli.command("gc-df-to-fasta")
@click.option(
    '--gc-hk-df',
    '--gc-df',
    type=click.Path(exists=True),
    required=True,
    help="Germinal center cell database - output \
            using the `wrangle_partis_annotation` command."
)
@click.option(
    '--header-col',
    '-h',
    multiple=True,
    type=str,
    default=["ID_HK"],
    help="A column from the germinal center df you \
            would like to concatinate to the header."
)
@click.option(
    '--sequence-col',
    '-s',
    multiple=True,
    type=str,
    default=["seq_nt_HC", "seq_nt_LC"],
    help="A column from the germinal center df you \
            would like to concatinate to the sequences."
)
@click.option(
    "--output",
    "-o",
    required=False,
    default="HK.fasta",
    help="click.Path to write the fasta.",
)
@click.option(
    '-n', '--add-naive',
    type=click.BOOL,
    default=True
)
def gc_df_to_fasta(gc_hk_df, header_col, sequence_col, output, add_naive):
    """
    A function to concatinate specified columns from
    the germinal center dataframe (output by
    `wrangle_partis_annotation` command) for both headers
    and sequences.
    """

    # gather the concatinated columns for headers
    gc_hk_df = pd.read_csv(gc_hk_df)
    headers = gc_hk_df[list(header_col)].apply(lambda x: "".join(x), axis=1)
    sequences = gc_hk_df[list(sequence_col)].apply(
        lambda x: "".join(x), axis=1)

    # write fasta
    with open(output, "w") as fasta:
        if add_naive:
            fasta.write(f">naive\n{naive_hk_bcr_nt}\n")
        for header, sequence in zip(headers, sequences):
            fasta.write(f">{header}\n{sequence}\n")

##### QUERY DATAFRAME
@cli.command("get-columns")
@click.option(
    '--dataframe',
    '-df',
    type=click.Path(exists=True),
    required=True,
    help="dataframe (csv) to query"
)
@click.option(
    '--columns',
    '-c',
    type=str,
    required=True,
    multiple=True,
    default=["ID_HK", "isotype_HC", "isotype_LC"],
    help='the key file for merging heavy and light chains'
)
@click.option(
    "--output",
    "-o",
    required=False,
    default="HK.fasta",
    help="click.Path to write the fasta.",
)
def get_columns(dataframe, columns, output):
    """
    Simply, a CLI wrapper for pandas DataFrame.loc on certain columns.
    """
    pd.read_csv(dataframe).loc[:, columns].to_csv(output, index=False)

##### FEATURIZE NODES
@cli.command("featurize-nodes")
@click.argument("gctree_file", type=click.Path(exists=True))
@click.argument("idmapfile", type=click.Path(exists=True))
@click.argument("variant_scores", type=click.Path(exists=False))
@click.option(
    "--naive_sites",
    type=click.Path(exists=False),
    default="https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
)
@click.argument("naive_sites", type=click.Path(exists=False))
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
@click.option("--tau_rei", type=float, default=0.5, help="Tau for REI computation.")
@click.option(
    "--output_dir",
    "-o",
    type=click.Path(exists=True),
    default=".",
    help="click.Path to write output files.",
)
@click.option(
    '-rid', '--render-idlabel',
    type=click.BOOL,
    default=True
)
def node_featurize(
    gctree_file,
    idmapfile,
    variant_scores,
    naive_sites,
    igk_idx,
    igh_frame,
    igk_frame,
    tau,
    tau0,
    tau_rei,
    output_dir,
    render_idlabel
):
    """Featurizes a gctree.CollapsedTree object using DMS data, outputing a new
    pickled tree object, a csv of node data, and rendered trees according to features.
    \b
    GCTREE_FILE: click.Path to pickled gctree.CollapsedTree object
    VARIANT_SCORES: click.Path to variant scores csv in DMS repository
    NAIVE_SITES: click.Path to sites csv in DMS repository
    """

    with open(idmapfile, "r") as fh:
        idmap = {}
        for line in fh:
            seqid, cell_ids = line.rstrip().split(",")
            idmap[seqid] = cell_ids.replace(":", ",")

    # DMS single mutant scores
    final_variant_scores = pd.read_csv(
        variant_scores, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
    )
    # remove linker sites
    final_variant_scores = final_variant_scores[final_variant_scores.chain != "link"]

    # position maps for scFv
    pos_df = pd.read_csv(
        naive_sites,
        dtype=dict(site=pd.Int16Dtype()),
        index_col="site_scFv",
    )
    # position maps for heavy and light chain that can be used in gctree render
    igh_pos_map = pos_df.loc[pos_df.chain
                             == "H", "site"].reset_index(drop=True)
    igk_pos_map = pos_df.loc[pos_df.chain
                             == "L", "site"].reset_index(drop=True)

    # load tree object
    with open(gctree_file, "rb") as f:
        tree = pickle.load(f)

    # generate REI node feature
    burst_stat(tree, "REI", tau=tau_rei)

    # generate LBI and LBR node features
    tree.local_branching(tau=tau, tau0=tau0)

    node_features = defaultdict(list)
    naive_igh_nt = tree.tree.sequence[:igk_idx]
    naive_igk_nt = tree.tree.sequence[igk_idx:]
    naive_igh_aa = aa(naive_igh_nt, igh_frame)
    naive_igk_aa = aa(naive_igk_nt, igk_frame)
    for node in tree.tree.traverse():
        igh_nt = node.sequence[:igk_idx]
        igk_nt = node.sequence[igk_idx:]
        igh_aa = aa(igh_nt, igh_frame)
        igk_aa = aa(igk_nt, igk_frame)
        igh_nt_mutations = nt_mutations(naive_igh_nt.upper(), igh_nt.upper())
        igk_nt_mutations = nt_mutations(naive_igk_nt.upper(), igk_nt.upper(), site_idx_offset=igk_idx)

        igh_mutations = mutations(naive_igh_aa, igh_aa, igh_pos_map, "(H)")
        igk_mutations = mutations(naive_igk_aa, igk_aa, igk_pos_map, "(L)")
        igh_has_stop = any("*" in x for x in igh_mutations)
        igk_has_stop = any("*" in x for x in igk_mutations)
        all_mutations = igh_mutations + igk_mutations
        node.add_feature("mutations", all_mutations)

        if len(node.isotype) == 1:
            isotype = next(iter(node.isotype.keys()))
        else:
            isotype = ",".join(
                f"{key}:{node.isotype[key]}" for key in node.isotype)

        node_features["name"].append(node.name)
        node_features["parent_name"].append(node.up.name if node.up else None)
        node_features["abundance"].append(node.abundance)
        node_features["sampled_cell_ids"].append(idmap[node.name] if node.name in idmap else "")

        node_features["n_nt_mutations_HC"].append(len(igh_nt_mutations))
        node_features["n_nt_mutations_LC"].append(len(igk_nt_mutations))
        node_features["IgH_nt_mutations"].append(" ".join(igh_nt_mutations))
        node_features["IgK_nt_mutations"].append(" ".join(igk_nt_mutations))

        node_features["n_aa_mutations_HC"].append(len(igh_mutations))
        node_features["n_aa_mutations_LC"].append(len(igk_mutations))
        node_features["IgH_aa_mutations"].append(" ".join(igh_mutations))
        node_features["IgK_aa_mutations"].append(" ".join(igk_mutations))

        node_features["IgH_productive"].append(not igh_has_stop)
        node_features["IgK_productive"].append(not igk_has_stop)
        node_features["isotype"].append(isotype)
        node_features["LBI"].append(node.LBI)
        node_features["LBR"].append(node.LBR)
        node_features["REI"].append(node.REI)
        node_features["descendant_abundance"].append(sum(descendant.abundance for descendant in node.traverse()))
        node_features["IgH_nt_sequence"].append(node.sequence[:igk_idx])
        node_features["IgH_aa_sequence"].append(str(igh_aa))
        node_features["IgK_nt_sequence"].append(node.sequence[igk_idx:])
        node_features["IgK_aa_sequence"].append(str(igk_aa))

        
        # substitutions to (hopefully) match those in FMVS
        ami = " ".join(all_mutations)
        node_features["aa_substitutions_IMGT"].append(ami)
        node.add_feature("aa_substitutions", ami)
        
        for i, phenotype in enumerate(["delta_bind_CGG", "delta_expr"]):

            # each of the additive scores from Tylers wrangling
            if final_variant_scores is not None:

                additive_score = ( 
                    np.nan if (igh_has_stop or igk_has_stop) 
                    else final_variant_scores.loc[all_mutations, phenotype].sum()
                )
                node_features[phenotype[:10]].append(additive_score)
                node.add_feature(phenotype[:10], additive_score)

        node.add_feature("delta_avidity", node.delta_bind + node.delta_expr)
        node_features["delta_avidity"].append(node.delta_avidity)

    df = pd.DataFrame(node_features).set_index("name")
    df.to_csv(f"{output_dir}/node_data.csv")

    tree.render(
        f"{output_dir}/mut_seq_annotated_nodes.svg",
        scale=20,
        idlabel=render_idlabel,
        frame=igh_frame,
        frame2=igk_frame,
        chain_split=igk_idx,
        position_map=igh_pos_map,
        position_map2=igk_pos_map,
    )

    # write the new featurized tree to a pickle file
    with open(f"{output_dir}/gctree.p", "wb") as f:
        pickle.dump(tree, f)


##### FEATURIZE SEQS
@cli.command("featurize-seqs")
@click.argument("hk_df", type=click.Path(exists=True))
@click.option("--igh_frame", type=int, default=1, help="Heavy chain reading frame.")
@click.option("--igk_frame", type=int, default=1, help="Light chain reading frame.")
@click.option(
    "--igk_idx",
    "-k",
    type=int,
    required=True,
    help="Start index of light chain in concatenated sequence.",
)
@click.option(
    "--variant_scores",
    type=click.Path(exists=False),
    default="https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/main/results/final_variant_scores/final_variant_scores.csv"
)
@click.option(
    "--naive_sites",
    type=click.Path(exists=False),
    default="https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
)
@click.option(
    "--output",
    "-o",
    type=click.Path(exists=False),
    default=".",
    help="click.Path to write output files.",
)
def featurize_seqs(
    hk_df,
    igh_frame,
    igk_frame,
    igk_idx,
    variant_scores,
    naive_sites,
    output
):
    """ Primarily for PB/MB cells - add DMS features to the dataframe
    including num mutations and additive Kd Scores for each sequence (row) in df

    HK_DF: A dataframe with paired HK seqs
    VARIANT_SCORES: click.Path to variant scores csv in DMS repository
    MULTI_VARIANT_SCORES: click.Path to final multi variant scores used for GE training
    NAIVE_SITES: click.Path to sites csv in DMS repository
    """
    hk_df = pd.read_csv(hk_df)

    # DMS single mutant scores, final from tyle
    final_variant_scores = pd.read_csv(
        variant_scores, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
    )
    # remove linker sites
    final_variant_scores = final_variant_scores[final_variant_scores.chain != "link"]
    
    # position maps for scFv
    pos_df = pd.read_csv(
        naive_sites,
        dtype=dict(site=pd.Int16Dtype()),
        index_col="site_scFv",
    )
    # position maps for heavy and light chain that can be used in gctree render
    igh_pos_map = pos_df.loc[pos_df.chain
                             == "H", "site"].reset_index(drop=True)
    igk_pos_map = pos_df.loc[pos_df.chain
                             == "L", "site"].reset_index(drop=True)


    # construct new dataframe with all sequence phenotype predictions
    seq_pheno_preds = defaultdict(list)
    naive_igh_nt = naive_hk_bcr_nt[:igk_idx]
    naive_igk_nt = naive_hk_bcr_nt[igk_idx:]

    naive_igh_aa = aa(naive_igh_nt, igh_frame)
    naive_igk_aa = aa(naive_igk_nt, igk_frame)
    
    # Well make a prediction for each of the observed sequences
    for idx, row in hk_df.iterrows():
        igh_nt = row.seq_nt_HC
        igk_nt = row.seq_nt_LC

        igh_nt_mutations = nt_mutations(naive_igh_nt.upper(), igh_nt.upper())
        igk_nt_mutations = nt_mutations(naive_igk_nt.upper(), igk_nt.upper(), site_idx_offset=igk_idx)

        # Infer mutations from nt seq available
        igh_aa = aa(row.seq_nt_HC, igh_frame)
        igk_aa = aa(row.seq_nt_LC, igk_frame)
        
        assert igh_aa == row.seq_aa_HC
        assert igk_aa == row.seq_aa_LC

        igh_mutations = mutations(naive_igh_aa, igh_aa, igh_pos_map, "(H)")
        igk_mutations = mutations(naive_igk_aa, igk_aa, igk_pos_map, "(L)")
        igh_has_stop = any("*" in x for x in igh_mutations)
        igk_has_stop = any("*" in x for x in igk_mutations)
        all_mutations = igh_mutations + igk_mutations


        # Now, we'll add all the phenotype predictions        
        # Individual chain mutations
        seq_pheno_preds["n_nt_mutations_HC"].append(len(igh_nt_mutations))
        seq_pheno_preds["n_nt_mutations_LC"].append(len(igk_nt_mutations))
        seq_pheno_preds["IgH_nt_mutations"].append(" ".join(igh_nt_mutations))
        seq_pheno_preds["IgK_nt_mutations"].append(" ".join(igk_nt_mutations))
        
        seq_pheno_preds["n_aa_mutations_HC"].append(len(igh_mutations))
        seq_pheno_preds["n_aa_mutations_LC"].append(len(igk_mutations))
        seq_pheno_preds["IgH_aa_mutations"].append(" ".join(igh_mutations))
        seq_pheno_preds["IgK_aa_mutations"].append(" ".join(igk_mutations))

        
        # substitutions to (hopefully) match those in FMVS
        seq_pheno_preds["aa_substitutions_IMGT"].append(" ".join(all_mutations))
        
        for i, phenotype in enumerate(["delta_bind_CGG", "delta_expr"]):

            if final_variant_scores is not None:

                seq_pheno_preds[phenotype[:10]].append(
                    np.nan if (igh_has_stop or igk_has_stop) 
                    else final_variant_scores.loc[all_mutations, phenotype].sum()
                )
                
        seq_pheno_preds["delta_avidity"].append(
            seq_pheno_preds["delta_bind"][-1] + seq_pheno_preds["delta_expr"][-1]
        )
        
    df = pd.DataFrame(seq_pheno_preds)
    ret = pd.concat([df, hk_df], axis=1)

    ret = ret.loc[:, [c for c in final_HK_col_order if c in ret.columns]]
    ret.to_csv(f"{output}", index=False)

if __name__ == '__main__':
    cli()
