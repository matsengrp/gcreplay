#!/usr/bin/env python
"""
@file: gcreplay-tools-dev.py

Right now a sandbox for throwing code into
it'll be a simple click CLI with all the logic right here.
"""
# TODO ADD all authors


import re
import pickle
import glob
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import torch
import click
from click import Choice, Path, command, group, option, argument
import regex

import matplotlib

from utils import *


#################################
# GLOBALS
#################################

# TODO we should probably just store the important sequences in fasta
# matter of fact, this already exists for the partis annotation. We'll
# just pass that.
naive_hk_bcr_nt = ("GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCT"
                   "CAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGA"
                   "TCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTA"
                   "CTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTAC"
                   "CTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATG"
                   "TCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAGACATTGTGATGACtCAGTCTCAAAAATT"
                   "CATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACT"
                   "AATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCT"
                   "ACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCAC"
                   "CATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCT"
                   "CTCACGTTCGGCTCGGGGACtAAGCTaGAAATAAAA")

naive_hk_bcr_aa = ("EVQLQESGPSLVKPSQTLSLTCSVTGDSITSGYWNWIRKFPGNKLEYMGYISYSG"
                   "STYYNPSLKSRISITRDTSKNQYYLQLNSVTTEDTATYYCARDFDVWGAGTTVTVSSDIVMTQS"
                   "QKFMSTSVGDRVSVTCKASQNVGTNVAWYQQKPGQSPKALIYSASYRYSGVPDRFTGSGSGTDF"
                   "TLTISNVQSEDLAEYFCQQYNSYPLTFGSGTKLEIK")

partis_airr_to_drop = [
    "j_sequence_end", "j_germline_start", "j_germline_end", "rev_comp",
    "junction", "clone_id", "vj_in_frame", "stop_codon", "np1", "np2", "duplicate_count",
    "cdr3_start", "cdr3_end", "cell_id", "v_support", "v_identity", "v_sequence_start",
    "v_sequence_end", "v_germline_start", "v_germline_end", "d_support", "d_identity",
    "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
    "j_support", "j_identity", "j_sequence_start", "j_sequence_end", "j_germline_start",
    "j_germline_end",
]

final_HK_col_order = [
    "ID_HK", "well", "HK_key_plate", "HK_key_mouse", "HK_key_gc", "HK_key_node", "HK_key_cell_type",
    "aa_substitutions_IMGT", "aa_sequence",

    "delta_bind_CGG_FVS_additive",
    "delta_expr_FVS_additive",
    "delta_psr_FVS_additive",

    "delta_bind_CGG_FMVS_additive",
    "delta_expr_FMVS_additive",
    "delta_psr_FMVS_additive",

    "delta_bind_CGG_FMVS_ground_truth",
    "delta_expr_FMVS_ground_truth",
    "delta_psr_FMVS_ground_truth",

    "delta_bind_CGG_tdms_model_pred",
    "delta_expr_tdms_model_pred",
    "delta_psr_tdms_model_pred",

    "delta_bind_CGG_tdms_linear_model_pred",
    "delta_expr_tdms_linear_model_pred",
    "delta_psr_tdms_linear_model_pred",

    "n_mutations_HC", "n_mutations_LC", "IgH_mutations", "IgK_mutations",
    "isotype_HC", "isotype_LC", "ID_HC", "ID_LC",
    "Productive_HC", "Productive_LC",
    "V_HC", "V_LC", "D_HC", "D_LC", "J_HC", "J_LC",
    "AAjunction_HC", "AAjunction_LC", "locus_HC", "locus_LC",
    "miseq_plate_HC", "miseq_plate_LC",
    "barcode_HC", "barcode_LC", "row_HC", "row_LC", "column_HC", "column_LC",
    "chain_HC", "chain_LC", "rank_HC", "rank_LC", "counts_HC", "counts_LC",
    "seq_nt_length_HC", "seq_nt_length_LC", "seq_aa_length_HC", "seq_aa_length_LC",
    "fasta_header_HC", "fasta_header_LC", "fasta_seq_HC", "fasta_seq_LC", 
    "partis_sequence_HC", "partis_sequence_LC",
    "seq_aa_HC", "seq_aa_LC",
    "seq_nt_HC", "seq_nt_LC", 
    # "date",
]


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


@cli.command("sort-fasta")
@click.option(
    '--fasta',
    type=Path(exists=True),
    required=True,
    help='a sequence count rank collapsed fasta with header format <rank>-<count>'
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Path to write the fasta.",
)
def sort_fasta(fasta: str, output: str):
    """ sort a fasta file by sorting the headers using
    pandas lexiconographical sort function """

    fasta_df = fasta_to_df(fasta).sort_values(by="id")
    with open(output, "w") as f2:
        for idx, row in fasta_df.iterrows():
            f2.write(f">{row.id}\n{row.seq}\n")


@cli.command("curate-high-count-seqs")
@click.option(
    '--fasta',
    type=Path(exists=True),
    required=True,
    help='a sequence count rank collapsed fasta with header format <rank>-<count>'
)
@click.option(
    '--count-threshold',
    '-ct',
    type=int,
    required=True,
    help='threshold of counts needed for a sequence to be included in output'
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Path to write the fasta.",
)
def curate_high_count_seqs(fasta: str, count_threshold: int, output: str):
    """Keep n sequences from  collapsed_fasta"""

    with open(fasta, "r") as f1, open(output, "w") as f2:
        for seq_record in SeqIO.parse(f1, 'fasta'):  # (generator)
            rank, count = str(seq_record.id).split("-")
            if int(count) < count_threshold: break
            f2.write(f">{seq_record.id}\n{seq_record.seq}\n")


@cli.command("wrangle-annotation")
@click.option(
    '--igh-airr',
    type=Path(exists=True),
    required=True,
    help='igh airr output from partis'
)
@click.option(
    '--igk-airr',
    type=Path(exists=True),
    required=True,
    help='igk airr output from partis'
)
@click.option(
    '--input-fasta',
    type=Path(exists=True),
    required=True,
    help='the fasta fed into partis for annotation'
)
@click.option(
    '--key-file',
    type=Path(exists=True),
    required=True,
    help='the key file for merging heavy and light chains'
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Path to write the fasta.",
)
def wrangle_annotation(
    igh_airr, igk_airr, input_fasta, key_file, output
):
    """
    Curate and wrangle the partis annotation of ranked BCR's
    for a single "PR" experiment.
    """

    # read in the airr formatted files into dataframes
    partis_igk = pd.read_csv(igk_airr, sep="\t")
    partis_igh = pd.read_csv(igh_airr, sep="\t")

    # filter out the failed annotations
    partis_igk.query("seqs_aa.notna()", engine="python", inplace=True)
    partis_igh.query("seqs_aa.notna()", engine="python", inplace=True)

    # TODO raise more useful error
    # explain that we expect some non-coding
    # nt at the end of each sequence
    for seq in partis_igh["seqs_aa"]:
        assert seq.endswith("X")

    # if the assertion is true we can simply
    # strip the 'X' from all aa sequences
    partis_igh.loc[:, "seq_aa"] = [s[:-1] for s in partis_igh["seqs_aa"]]
    partis_igh = partis_igh.drop("seqs_aa", axis=1)
    partis_igk = partis_igk.rename({"seqs_aa": "seq_aa"}, axis=1)

    # merge the igk and igh annotations
    # partis_airr = partis_igk.append(partis_igh)
    partis_airr = pd.concat([partis_igk, partis_igh])
    partis_airr.loc[:, "seq_nt"] = [seq.lower()
                                    for seq in partis_airr["sequence"]]

    # drop unnecessary columns
    partis_airr.drop(partis_airr_to_drop, axis=1, inplace=True)


    # parse the fasta to merge in the important information
    # TODO add plate int to this as 'barcode', keep plate for the hell of it
    parsed_input_fasta = bcr_fasta_to_df(input_fasta, parse_nextflow_header)

    # Merge in the parsed fasta columns
    # TODO once we're getting rid of unmatched we can do a left join? assert this.
    partis_airr = partis_airr.merge(
        parsed_input_fasta, how="inner", on="sequence_id")

    partis_airr["ID"] = [
        f"{r.ngs_date}{r.plate}{r.well}{r.chain}"
        for i, r in partis_airr.iterrows()
    ]

    partis_airr = partis_airr.rename(
        {
            "v_call": "V",
            "d_call": "D",
            "j_call": "J",
            "productive": "Productive",
            "junction_aa": "AAjunction",
            "sequence_id": "fasta_header",
            "plate" : "miseq_plate",
            "sequence" : "partis_sequence"
        },
        axis=1
    )

    # compute lengths
    partis_airr.loc[:, "seq_nt_length"] = [
        len(seq) for seq in partis_airr["seq_nt"]]
    partis_airr.loc[:, "seq_aa_length"] = [
        len(seq) for seq in partis_airr["seq_aa"]]

    query_string = "(locus == 'IGH' & seq_nt_length == 337) | (locus == 'IGK' & seq_nt_length == 321)"

    # Trim the rest of partis IGH
    partis_airr.query(query_string, engine="python", inplace=True)
    seq_nts = []
    for idx, row in partis_airr.iterrows():
        seq_nts.append(row.seq_nt[:-1] if row.locus == "IGH" else row.seq_nt)
    partis_airr.loc[:, "seq_nt"] = seq_nts

    # re-compute lengths
    partis_airr.loc[:, "seq_nt_length"] = [
        len(seq) for seq in partis_airr["seq_nt"]]
    partis_airr.loc[:, "seq_aa_length"] = [
        len(seq) for seq in partis_airr["seq_aa"]]

    partis_airr.loc[:, "isotype"] = [
        infer_igh_isotypes(row.fasta_seq, 3)
        if row.locus == "IGH" else "IgK"
        for idx, row in partis_airr.iterrows()
    ]
    no_iso = partis_airr.loc[partis_airr.isotype.isnull()].index
    partis_airr.drop(no_iso, inplace=True)

    # Some cells are spread across wells.
    # This means we need to merge the counts across plates 
    # where wells are identical

    # TODO This is kinda dirty
    # I feel like going forward we should enforce a single cell per well

    print(partis_airr)
    assert len(partis_airr) == len(set(partis_airr.index))

    key_file_df = pd.read_csv(key_file)
    for idx, row in key_file_df.iterrows():
        for chain, column in zip(["IGH" ,"IGK"], ["hc_barcode", "lc_barcode"]):
            key_chain_barcodes = [int(c) for c in str(row[column]).split(".")]

            if len(key_chain_barcodes) > 1:
        
                # Let's make sure this code isn't effecting any unwanted
                # normal GC's. I assume it won't
                assert row.gc == 20

                key_row = row.row.split(".")
                key_col = [int(c) for c in row.col.split(".")]
                q = (
                    f"(locus == '{chain}') & "
                    f"(barcode.isin({key_chain_barcodes})) & "
                    f"(column.isin({key_col})) & "
                    f"(row.isin({key_row}))"
                )

                C = partis_airr.query(q, engine="python")

                to_throw = set(C.index.values)
                for well, well_multiplate_df in C.groupby("well"):

                    well_summed_seqs = well_multiplate_df.groupby("fasta_seq").sum()
                    sorted_sum_well_seqs = well_summed_seqs.sort_values("counts", ascending=False)
                    winning_seq = sorted_sum_well_seqs.index.values[0]
                    winning_counts = sorted_sum_well_seqs.counts.values[0]
                    winning_original_entries = well_multiplate_df.query(f"fasta_seq == '{winning_seq}'")
                    winning_index = winning_original_entries.index.values[0]
                    to_throw.remove(winning_index)
                    partis_airr.loc[winning_index, "counts"] = winning_counts
                    partis_airr.loc[winning_index, "rank"] = 1

                # now drop all the non winning BCR's
                partis_airr.drop(to_throw, axis=0, inplace=True)

    # TODO shall we export those which are not rank 1, count 10?
    # TODO make these parameters. I guess rank one is kind of innevitable
    partis_airr.query("(rank == 1)", engine="python", inplace=True)
    for well, well_df in partis_airr.groupby("ID"):
        assert len(well_df) == 1

    # Merge heavy and light chain
    GC_df = merge_heavy_light_chains(partis_airr, key_file_df)

    # TODO This should be somewhere else, but it's a qc on number of Ns in a sequence.
    throw_seq = [idx for idx, row in GC_df.iterrows() if (
        row.seq_nt_HC + row.seq_nt_LC).count("n") > 30]
    GC_df.drop(throw_seq, axis=0, inplace=True)

    # clean up
    # TODO date?
    ret = GC_df.loc[:,[c for c in final_HK_col_order if c in GC_df.columns]]
    ret.to_csv(output, index=False)


##### QUERY DATAFRAME
# TODO should the query be done here and kept in memory to
# save the size of intermediate files?
@cli.command("gc-df-to-fasta")
@click.option(
    '--gc-hk-df',
    '--gc-df',
    type=Path(exists=True),
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
    help="Path to write the fasta.",
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


##### SAMPLE DATAFRAME
@cli.command("sample-10-df")
@click.option(
    '--dataframe',
    '-df',
    type=Path(exists=True),
    required=True,
    help="dataframe (csv) to query"
)
@click.option(
    "--output",
    "-o",
    required=True,
)
def sample_10_df(dataframe, output):
    """
    Simply, a CLI wrapper for pandas DataFrame.sample on certain columns.
    """
    df = pd.read_csv(dataframe)
    df.sample(min(10, len(df))).to_csv(output, index=False)



##### GROUBY DATAFRAME
@cli.command("df-groupby")
@click.option(
    '--dataframe',
    '-df',
    type=Path(exists=True),
    required=True,
    help="dataframe (csv) to query"
)
@click.option(
    '--columns',
    '-c',
    type=str,
    required=True,
    multiple=True,
    default=["HK_key_mouse", "HK_key_node", "HK_key_gc", "HK_key_cell_type"],
    help='the key file for merging heavy and light chains'
)
@click.option("--sample", type=int, default=None, help="Heavy chain reading frame.")
@click.option(
    "--output-prefix",
    "-o",
    required=False,
    default="grouped",
    help="prefix of filename before group string",
)
def df_groupby(dataframe, columns, sample, output_prefix):
    """
    Simply, a CLI wrapper for pandas DataFrame.groupby on certain columns.
    """
    df = pd.read_csv(dataframe)
    for group, groupdf in df.groupby(list(columns)):
        group_string = "-".join([str(gi) for gi in group])
        if sample != None:
            groupdf.sort_values(by="ID_HK").sample(min(10,len(groupdf)), random_state=1).to_csv(f"{output_prefix}-{group_string}.csv", index=False)
        else:
            groupdf.sort_values(by="ID_HK").to_csv(f"{output_prefix}-{group_string}.csv", index=False)


##### QUERY DATAFRAME
@cli.command("get-columns")
@click.option(
    '--dataframe',
    '-df',
    type=Path(exists=True),
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
    help="Path to write the fasta.",
)
def get_columns(dataframe, columns, output):
    """
    Simply, a CLI wrapper for pandas DataFrame.loc on certain columns.
    """
    pd.read_csv(dataframe).loc[:, columns].to_csv(output, index=False)


##### QUERY DATAFRAME
@cli.command("query-df")
@click.option(
    '--dataframe',
    '-df',
    type=Path(exists=True),
    required=True,
    help="dataframe (csv) to query"
)
@click.option(
    '--query-string',
    '-q',
    type=str,
    required=True,
    help='the key file for merging heavy and light chains'
)
@click.option(
    "--output",
    "-o",
    required=False,
    default="HK.fasta",
    help="Path to write the fasta.",
)
def query_df(dataframe, query_string, output):
    """
    Simply, a CLI wrapper for pandas DataFrame.query.
    """
    pd.read_csv(dataframe).query(query_string).to_csv(output, index=False)


##### FEATURIZE NODES
@cli.command("featurize-nodes")
@click.argument("gctree_file", type=click.Path(exists=True))
@click.argument("idmapfile", type=click.Path(exists=True))
@click.argument("variant_scores", type=click.Path(exists=False))
@click.option(
    "--multi_variant_scores",
    type=click.Path(exists=False),
    default=None
)
@click.option(
    "--tdms_model",
    type=click.Path(exists=False),
    default=None
)
@click.option(
    "--tdms_model_linear",
    type=click.Path(exists=False),
    default=None
)
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
@click.option(
    "--output_dir",
    "-o",
    type=click.Path(exists=True),
    default=".",
    help="Path to write output files.",
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
    multi_variant_scores,
    tdms_model,
    tdms_model_linear,
    naive_sites,
    igk_idx,
    igh_frame,
    igk_frame,
    tau,
    tau0,
    output_dir,
    render_idlabel
):
    """Featurizes a gctree.CollapsedTree object using DMS data, outputing a new
    pickled tree object, a csv of node data, and rendered trees according to features.
    \b
    GCTREE_FILE: Path to pickled gctree.CollapsedTree object
    VARIANT_SCORES: Path to variant scores csv in DMS repository
    NAIVE_SITES: Path to sites csv in DMS repository
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

    # Final Multi Variant Scores 
    # TODO Pass this in pipeline
    final_multi_variant_scores = (
        pd.read_csv(multi_variant_scores, index_col="aa_substitutions_IMGT")
        #.query("library == 'dms'") # TODO we are not going to add pred with octet data
        if multi_variant_scores is not None
        else None
    )

    # torchdms models
    tdms_model_binary = (
        torch.load(tdms_model, map_location="cpu")
        if tdms_model is not None
        else None
    )

    tdms_model_linear_binary = (
        torch.load(tdms_model_linear, map_location="cpu")
        if tdms_model_linear is not None
        else None
    )

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

    # generate LBI and LBR node features
    tree.local_branching(tau=tau, tau0=tau0)

    node_features = defaultdict(list)
    naive_igh_aa = aa(tree.tree.sequence[:igk_idx], igh_frame)
    naive_igk_aa = aa(tree.tree.sequence[igk_idx:], igk_frame)
    for node in tree.tree.traverse():
        igh_aa = aa(node.sequence[:igk_idx], igh_frame)
        igk_aa = aa(node.sequence[igk_idx:], igk_frame)

        # Make the aa seq for tdms 
        aa_tdms = pos_df.amino_acid.copy()
        aa_tdms.iloc[pos_df.chain == "H"] = igh_aa
        # note: replay light chains are shorter than dms seq by one aa
        aa_tdms.iloc[(pos_df.chain == "L") & (pos_df.index < pos_df.index[-1])] = igk_aa
        aa_tdms_seq = "".join(aa_tdms)
        node.add_feature("aa_sequence", aa_tdms_seq)

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
        node_features["n_mutations_HC"].append(len(igh_mutations))
        node_features["n_mutations_LC"].append(len(igk_mutations))
        node_features["IgH_mutations"].append(",".join(igh_mutations))
        node_features["IgK_mutations"].append(",".join(igk_mutations))
        node_features["IgH_productive"].append(not igh_has_stop)
        node_features["IgK_productive"].append(not igk_has_stop)
        node_features["isotype"].append(isotype)
        node_features["LBI"].append(node.LBI)
        node_features["LBR"].append(node.LBR)
        node_features["descendant_abundance"].append(sum(descendant.abundance for descendant in node.traverse()))
        node_features["IgH_nt_sequence"].append(node.sequence[:igk_idx])
        node_features["IgH_aa_sequence"].append(str(igh_aa))
        node_features["IgK_nt_sequence"].append(node.sequence[igk_idx:])
        node_features["IgK_aa_sequence"].append(str(igk_aa))
        node_features["aa_sequence"].append(aa_tdms_seq)

        
        # substitutions to (hopefully) match those in FMVS
        ami = " ".join(all_mutations)
        node_features["aa_substitutions_IMGT"].append(ami)
        node.add_feature("aa_substitutions", ami)
        
        for i, phenotype in enumerate(["delta_bind_CGG", "delta_expr", "delta_psr"]):

            # each of the additive scores from Tylers wrangling
            if final_variant_scores is not None:

                additive_score = ( 
                    np.nan if (igh_has_stop or igk_has_stop) 
                    else final_variant_scores.loc[all_mutations, phenotype].sum()
                )
                node_features[f"{phenotype}_FVS_additive"].append(additive_score)
                node.add_feature(f"{phenotype}_FVS_additive", additive_score)


            # each of the 'additive' phenotypes from GE training prep
            if final_multi_variant_scores is not None:
                
                additive_score = (
                    np.nan if (igh_has_stop or igk_has_stop)
                    else final_multi_variant_scores.loc[all_mutations, phenotype].sum()
                )

                node_features[f"{phenotype}_FMVS_additive"].append(additive_score)
                node.add_feature(f"{phenotype}_FMVS_additive", additive_score)

            # each of the phenotype predictions for a global epistasis non linear model
            if tdms_model_binary is not None:
                seq_binary = tdms_model_binary.seq_to_binary(aa_tdms_seq)
                pred = tdms_model_binary(seq_binary).detach().numpy()[i]
                node_features[f"{phenotype}_tdms_model_pred"].append(pred)
                node.add_feature(f"{phenotype}_tdms_model_pred", pred)

                # TODO call tdms_model.to_latent ... only if the model latent makes sense
                # for each of the latent nodes attached to a prediction individually.
                #seq_binary = tdms_model_binary.seq_to_binary(aa_tdms_seq)
                #node_features[f"{phenotype}_tdms_model_pred_latent"] = (
                #    tdms_model_binary.to_latent(seq_binary).detach().numpy()[i]
                #)

            # each of the predictions from a perceptron tdms model
            if tdms_model_linear_binary is not None:
                seq_binary = tdms_model_linear_binary.seq_to_binary(aa_tdms_seq)
                pred = tdms_model_linear_binary(seq_binary).detach().numpy()[i]
                node_features[f"{phenotype}_tdms_linear_model_pred"].append(pred)
                node.add_feature(f"{phenotype}_tdms_linear_model_pred", pred)

    df = pd.DataFrame(node_features).set_index("name")
    df.to_csv(f"{output_dir}/node_data.csv")

    phenotypes = [
        "delta_bind_CGG_FVS_additive", 
        "delta_expr_FVS_additive", 
        "delta_psr_FVS_additive"
    ]
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
        colormap = tree.feature_colormap(
            phenotype, cmap=cmap, vmin=vmin, vmax=vmax)
        # need this loop to render to svg and notebook
        tree.render(
            f"{output_dir}/{phenotype}.svg",
            # scale=None, branch_margin=-7,
            scale=20,
            idlabel=render_idlabel,
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
    "--multi_variant_scores",
    type=click.Path(exists=False),
    default=None
)
@click.option(
    "--tdms_model",
    type=click.Path(exists=False),
    default=None
)
@click.option(
    "--tdms_model_linear",
    type=click.Path(exists=False),
    default=None
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
    help="Path to write output files.",
)
def featurize_seqs(
    hk_df,
    igh_frame,
    igk_frame,
    igk_idx,
    variant_scores,
    multi_variant_scores,
    tdms_model,
    tdms_model_linear,
    naive_sites,
    output
):
    """ Primarily for PB/MB cells - add DMS features to the dataframe
    including num mutations and additive Kd Scores for each sequence (row) in df

    HK_DF: A dataframe with paired HK seqs
    VARIANT_SCORES: Path to variant scores csv in DMS repository
    MULTI_VARIANT_SCORES: Path to final multi variant scores used for GE training
    NAIVE_SITES: Path to sites csv in DMS repository
    """
    hk_df = pd.read_csv(hk_df)

    # DMS single mutant scores, final from tyle
    final_variant_scores = pd.read_csv(
        variant_scores, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
    )
    # remove linker sites
    final_variant_scores = final_variant_scores[final_variant_scores.chain != "link"]

    # Final Multi Variant Scores 
    final_multi_variant_scores = (
        pd.read_csv(multi_variant_scores, index_col="aa_substitutions_IMGT")
        # .query("library == 'dms'") # TODO we are not going to add pred with octet data
        if multi_variant_scores is not None
        else None
    )

    # torchdms models
    tdms_model_binary = (
        torch.load(tdms_model, map_location="cpu")
        if tdms_model is not None
        else None
    )

    tdms_model_linear_binary = (
        torch.load(tdms_model_linear, map_location="cpu")
        if tdms_model_linear is not None
        else None
    )

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
    naive_igh_aa = aa(naive_hk_bcr_nt[:igk_idx], igh_frame)
    naive_igk_aa = aa(naive_hk_bcr_nt[igk_idx:], igk_frame)

    # Well make a prediction for each of the observed sequences
    for idx, row in hk_df.iterrows():

        # Infer mutations from nt seq available
        igh_aa = aa(row.seq_nt_HC, igh_frame)
        igk_aa = aa(row.seq_nt_LC, igk_frame)

        # Make the aa seq for tdms 
        aa_tdms = pos_df.amino_acid.copy()
        aa_tdms.iloc[pos_df.chain == "H"] = igh_aa
        # note: replay light chains are shorter than dms seq by one aa
        aa_tdms.iloc[(pos_df.chain == "L") & (pos_df.index < pos_df.index[-1])] = igk_aa
        aa_tdms_seq = "".join(aa_tdms)
        seq_pheno_preds["aa_sequence"].append(aa_tdms_seq)

        igh_mutations = mutations(naive_igh_aa, igh_aa, igh_pos_map, "(H)")
        igk_mutations = mutations(naive_igk_aa, igk_aa, igk_pos_map, "(L)")
        igh_has_stop = any("*" in x for x in igh_mutations)
        igk_has_stop = any("*" in x for x in igk_mutations)
        all_mutations = igh_mutations + igk_mutations


        # Now, we'll add all the phenotype predictions        
        # Individual chain mutations
        seq_pheno_preds["IgH_mutations"].append(igh_mutations)
        seq_pheno_preds["IgK_mutations"].append(igk_mutations)
        
        # substitutions to (hopefully) match those in FMVS
        seq_pheno_preds["aa_substitutions_IMGT"].append(" ".join(all_mutations))
        
        for i, phenotype in enumerate(["delta_bind_CGG", "delta_expr", "delta_psr"]):

            if final_variant_scores is not None:

                seq_pheno_preds[f"{phenotype}_FVS_additive"].append(
                    np.nan if (igh_has_stop or igk_has_stop) 
                    else final_variant_scores.loc[all_mutations, phenotype].sum()
                )

            # each of the 'additive' phenotypes from GE training prep
            if final_multi_variant_scores is not None:

                # compute additive score from final multi variant scores
                additive_score = (
                    np.nan if (igh_has_stop or igk_has_stop) 
                    else final_multi_variant_scores.loc[all_mutations, phenotype].sum()
                )

                seq_pheno_preds[f"{phenotype}_FMVS_additive"].append(additive_score)

            # each of the phenotype predictions for a global epistasis non linear model
            if tdms_model_binary is not None:
                seq_binary = tdms_model_binary.seq_to_binary(aa_tdms_seq)
                seq_pheno_preds[f"{phenotype}_tdms_model_pred"].append(
                    tdms_model_binary(seq_binary).detach().numpy()[i]
                )

                # TODO call tdms_model.to_latent ... only if the model latent makes sense
                # for each of the latent nodes attached to a prediction individually.
                #seq_binary = tdms_model_binary.seq_to_binary(aa_tdms_seq)
                #seq_pheno_preds[f"{phenotype}_tdms_model_pred_latent"] = (
                #    tdms_model_binary.to_latent(seq_binary).detach().numpy()[i]
                #)

            # each of the predictions from a perceptron tdms model
            if tdms_model_linear_binary is not None:
                seq_binary = tdms_model_linear_binary.seq_to_binary(aa_tdms_seq)
                seq_pheno_preds[f"{phenotype}_tdms_linear_model_pred"].append(
                    tdms_model_linear_binary(seq_binary).detach().numpy()[i]
                )
        
    df = pd.DataFrame(seq_pheno_preds)
    ret = pd.concat([df, hk_df], axis=1)

    # TODO ADD NEW PHENOTYPES
    ret = ret.loc[:, [c for c in final_HK_col_order if c in ret.columns]]
    ret.to_csv(f"{output}", index=False)



##### MERGE FEATURIZED TREE SEQS
@cli.command("merge-results")
@click.option("--glob-pattern", type=str, default="PR*")
@click.option(
    "--output-gcdf",
    "-ogcdf",
    type=click.Path(exists=False),
    default="observed-seqs.csv",
    help="Path to write output files.",
)
@click.option(
    "--output-gctree",
    "-ogctree",
    type=click.Path(exists=False),
    default="gctree-node-data.csv",
    help="Path to write output files.",
)
def merge_results(
    glob_pattern,
    output_gcdf,
    output_gctree
):
    """
    merge all the gc results.
    """
    
    gcdf, gctree = pd.DataFrame(), pd.DataFrame()
    for gct_out in glob.glob(glob_pattern):
        gcdf = pd.concat([gcdf, pd.read_csv(f"{gct_out}/observed_seqs.csv")])
        PR, mouse, node, gc, ct = gct_out.split("-")

        if ct == "GC":
            node_data = pd.read_csv(f"{gct_out}/node_data.csv")
            node_data.insert(0, "PR", PR)
            node_data.insert(1, "HK_key_mouse", mouse)
            node_data.insert(2, "HK_key_node", node)
            node_data.insert(3, "HK_key_gc", gc)
            node_data.insert(4, "HK_key_cell_type", ct)
            gctree = pd.concat([gctree, node_data])

    gcdf.to_csv(output_gcdf, index=False)
    gctree.to_csv(output_gctree, index=False)


# TODO We need this until we get setup.py
# for a real mf python package.
if __name__ == '__main__':
    cli()
