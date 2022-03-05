#!/usr/bin/env python
"""
@file: gcreplay-tools

Right now a sandbox for throwing code into
it'll be a simple click CLI with all the logic right here.
"""
# TODO ADD all authors


import re
import pickle

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import click
from click import Choice, Path, command, group, option, argument
import regex

import matplotlib
#matplotlib.use("Qt5Agg")


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


#################################
# HELPERS
#################################

def bcr_fasta_to_df(fasta_fp, id_parse_fn, **kwargs):
    """convert a fasta file pointer to dataframe after
    parsing the id with some function returning
    the columns defined (less the sequence column)"""

    columns = [
        'sequence_id',
        'plate',
        'barcode',
        'well',
        'row',
        'column',
        'chain',
        'rank',
        'counts',
        'seq_input'
    ]

    ret = pd.DataFrame({c: [] for c in columns})
    with open(fasta_fp) as fasta_file:  # Will close handle cleanly
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            bcr_meta = id_parse_fn(seq_record.id)
            if bcr_meta != -1:
                bcr_meta["seq_input"] = str(seq_record.seq)
                ret = ret.append(pd.Series(bcr_meta), ignore_index=True)
    return ret


def parse_nextflow_header(header: str):
    """parse fasta header and return rank, counts, well, plate, and chain,
    we do not expect the `>` to be included in the header"""
    if "unmatched" in header:
        return -1

    pr, date, plate, well, chain, _, rank_count = header.split(".")
    bcr_ranking, bcr_count = rank_count.split("-")

    return {
        "sequence_id": header,
        "date": date,
        "plate": plate,
        "barcode": int(plate[1:]),
        "well": well,
        "row": well[0],
        "column": int(well[1:]),
        "chain": chain,
        "rank": int(bcr_ranking),
        "counts": int(bcr_count)
    }


def plot_venn_stub(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    feature_groups: list,
    out="venn.png",
    **kwargs
):
    """
    plot a venn diagram of differences between grouped column features
    of two dataframes
    """

    #fig, ax = plt.subplots(figsize=[6,6])
    #unique_bcr_qualities = ["plate", "well","chain", "sequence"]
    #nf_well_bcrs = set(nextflow_input_bcrs.groupby(unique_bcr_qualities).groups)
    #tat_well_bcrs = set(tatsuya_input_bcrs.groupby(unique_bcr_qualities).groups)
    #v = venn2([nf_well_bcrs, tat_well_bcrs], ["nextflow", "tatsuya"], ax=ax)
    #ax.set_title("venn diagram of well/plate/chain specific \nranked BCR seqs from both methods\npre-annotation")
    #fig.savefig("input-sequence-differences.png")
    #plt.tight_layout()
    #plt.show()

    pass


def test_sequence_counts_stub():
    """
    TODO
    """

    #in_both = set.intersection(nf_well_bcrs, tat_well_bcrs)
    #for rkd_bcr in iter(in_both):
    #    queries = []
    #    for i, attr in enumerate(unique_bcr_qualities):
    #        queries.append(f"({attr} == '{rkd_bcr[i]}')")
    #    query_string = " & ".join(queries)
    #    nf_bcr = nextflow_input_bcrs.query(query_string, engine='python')
    #    tat_bcr = tatsuya_input_bcrs.query(query_string, engine='python')
    #    assert nf_bcr["count"].values[0] == tat_bcr["count"].values[0]
    #print("success! Sequences that appear in both methods have the exact same count!")

    pass


# TODO
def infer_igh_isotypes(bcr_sequence, num_mm=1):
    """Use fuzzy string matching to infer the isotype of bcr sequences
    in the cell database - allowing for `num_mm` mismatches. This function
    will modify the passed in dataframe and will return the first match it finds
    """

    # Motifs associated with full input BCR including constant region
    motif_map = {
        "IgM": ["atgtcttccccct"],
        "IgG1": ["atggtgaccctggg"],
        "IgG2": ["ggctcctcggtgactcta", "tcggtgactctagg", "gtgtggagatacaactgg"],
        "IgG3": ["cccttggtccctggctgcggtgacacat", "cacatctggatcctcggtgaca"],
        "IgA": ["tctgcgagaaatcccaccatcta"]
    }

    # match and return first valid isotype
    for key, value in motif_map.items():
        for motif in value:
            m = regex.findall("({m}){{e<={nm}}}".format(
                m=motif, nm=num_mm), bcr_sequence.lower())
            if len(m) > 0:
                return key
    # TODO Log
    print(f"WARNING: no isotype match for: {bcr_sequence}")

    # if no istype motif is found
    return np.nan


def _test_infer_igh_isotypes() -> None:
    assert infer_igh_isotypes("cctgaagtctgcgagaaatcccaccatctatctta") == "IgA"
    assert infer_igh_isotypes("cctgaagtctgcgagatatcccaccatctatctta") == "IgA"


def aa(sequence, frame):
    """Amino acid translation of nucleotide sequence in frame 1, 2, or 3."""
    return Seq(
        sequence[(frame - 1): (frame - 1
                               + (3 * ((len(sequence) - (frame - 1)) // 3)))]
    ).translate()


def mutations(naive_aa, aa, pos_map, chain_annotation):
    """Amino acid substitutions between two sequences, in IMGT coordinates."""
    return [
        f"{aa1}{pos_map[pos]}{chain_annotation}{aa2}"
        for pos, (aa1, aa2) in enumerate(zip(naive_aa, aa))
        if aa1 != aa2
    ]


def merge_heavy_light_chains(
        cell_df: pd.DataFrame,
        key_file: pd.DataFrame
) -> pd.DataFrame:
    """
    """

    GC_df = pd.DataFrame()
    for idx, row in key_file.iterrows():
        key_row = row.row.split(".")
        key_col = [int(c) for c in row.col.split(".")]

        HC = cell_df.query(
            f"(locus == 'IGH') & (barcode == {row.hc_barcode}) & (column.isin({key_col})) & (row.isin({key_row}))",
            engine="python"
        )
        LC = cell_df.query(
            f"(locus == 'IGK') & (barcode == {row.lc_barcode}) & (column.isin({key_col})) & (row.isin({key_row}))",
            engine="python"
        )

        GC = HC.merge(LC, on="well", suffixes=("_HC", "_LC"))
        GC["mouse_HC"] = row.mouse
        GC["GC_num_HC"] = row.gc
        GC["node_HC"] = row.node
        GC["cell_type_HC"] = row.cell_type
        GC["plate_num_HC"] = row.plate

        GC_df = pd.concat([GC_df, GC])

    # TODO we're missing just 2 mismatches compared to tatsuya, I believe - check this
    GC_df.loc[:, "ID_HK"] = [f"{i}K" for i in GC_df["ID_HC"]]
    GC_df = GC_df.query("chain_HC != chain_LC")

    return GC_df


#################################
# CLI
#################################

"""
# Usage example
python gcreplay-tools.py wrangle-annotation \
        --igh-airr 2022-02-09/partis_annotation/PR-1-6/engrd/single-chain/partition-igh.tsv \
        --igk-airr 2022-02-09/partis_annotation/PR-1-6/engrd/single-chain/partition-igk.tsv \
        --input-fasta 2022-02-09/ranked_bcr_sequences_per_well/PR-1-6.fasta
"""


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

# TODO option, filter unproductive cells
# TODO option, mm in isotype inference
# TODO option, counts threshold
# entry point


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
    partis_airr = partis_igk.append(partis_igh)
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
        f"{r.date}{r.plate}{r.well}{r.chain}"
        for i, r in partis_airr.iterrows()
    ]

    partis_airr = partis_airr.rename(
        {
            "v_call": "V",
            "d_call": "D",
            "j_call": "J",
            "productive": "Productive",
            "junction_aa": "AAjunction",
            "sequence_id": "fasta_header"
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
        infer_igh_isotypes(row.seq_input)
        if row.locus == "IGH" else "IgK"
        for idx, row in partis_airr.iterrows()
    ]

    # TODO shall we export those which are not rank 1, count 10?
    # TODO make these parameters. I guess rank one is kind of innevitable
    partis_airr.query("(rank == 1) & (counts >= 10)",
                      engine="python", inplace=True)
    for well, well_df in partis_airr.groupby("ID"):
        assert len(well_df) == 1

    GC_df = merge_heavy_light_chains(partis_airr, pd.read_csv(key_file))

    # TODO This should be somewhere else, but it's a qc on number of Ns in a sequence.
    throw_seq = [idx for idx, row in GC_df.iterrows() if (
        row.seq_nt_HC + row.seq_nt_LC).count("n") > 30]
    GC_df.drop(throw_seq, axis=0, inplace=True)

    GC_df.to_csv(output, index=False)


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
    default=["mouse_HC", "GC_num_HC", "cell_type_HC"],
    help='the key file for merging heavy and light chains'
)
@click.option(
    "--output-prefix",
    "-o",
    required=False,
    default="grouped",
    help="prefix of filename before group string",
)
def df_groupby(dataframe, columns, output_prefix):
    """
    Simply, a CLI wrapper for pandas DataFrame.groupby on certain columns.
    """
    df = pd.read_csv(dataframe)
    for group, groupdf in df.groupby(list(columns)):
        group_string = "-".join([str(gi) for gi in group])
        groupdf.to_csv(f"{output_prefix}-{group_string}.csv", index=False)


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
    '--max-n',
    '-n',
    type=str,
    required=True,
    help='max number of ns allowed in seq_nt_HC'
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
    output_dir,
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
    igh_pos_map = pos_df.loc[pos_df.chain
                             == "H", "site"].reset_index(drop=True)
    igk_pos_map = pos_df.loc[pos_df.chain
                             == "L", "site"].reset_index(drop=True)

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
        igh_has_stop = any("*" in x for x in igh_mutations)
        igk_has_stop = any("*" in x for x in igk_mutations)
        all_mutations = igh_mutations + igk_mutations
        node.add_feature("mutations", all_mutations)

        if len(node.isotype) == 1:
            isotype = next(iter(node.isotype.keys()))
        else:
            isotype = ",".join(
                f"{key}:{node.isotype[key]}" for key in node.isotype)
        row = [
            node.name,
            node.up.name if node.up else None,
            node.abundance,
            idmap[node.name] if node.name in idmap else "",
            node.sequence[:igk_idx],
            str(igh_aa),
            node.sequence[igk_idx:],
            str(igk_aa),
            ",".join(igh_mutations),
            ",".join(igk_mutations),
            not igh_has_stop,
            not igk_has_stop,
            isotype,
            node.LBI,
            node.LBR,
        ]
        for phenotype in phenotypes:
            node.add_feature(
                phenotype,
                np.nan if (igh_has_stop or igk_has_stop) else dms_df.loc[all_mutations, phenotype].sum(
                ),
            )
            row.append(getattr(node, phenotype))
        dat.append(row)

    # write dataframe
    columns = [
        "name",
        "parent_name",
        "abundance",
        "sampled_cell_ids",
        "IgH_nt_sequence",
        "IgH_aa_sequence",
        "IgK_nt_sequence",
        "IgK_aa_sequence",
        "IgH_mutations",
        "IgK_mutations",
        "IgH_productive",
        "IgK_productive",
        "isotype",
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
        colormap = tree.feature_colormap(
            phenotype, cmap=cmap, vmin=vmin, vmax=vmax)
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


# WSD note: pending discussion, but I think we may not need this, so commenting out
##### FEATURIZE SEQS
# @cli.command("featurize-seqs")
# @click.argument("hk_df", type=click.Path(exists=True))
# @click.argument(
#     "variant_scores",
#     type=click.Path(exists=False),
#     default="https://media.githubusercontent.com/media/jbloomlab/Ab-CGGnaive_DMS/main/results/final_variant_scores/final_variant_scores.csv"
# )
# @click.argument(
#     "naive_sites",
#     type=click.Path(exists=False),
#     default="https://raw.githubusercontent.com/jbloomlab/Ab-CGGnaive_DMS/main/data/CGGnaive_sites.csv"
# )
# @click.option(
#     "--output",
#     "-o",
#     type=click.Path(exists=True),
#     default=".",
#     help="Path to write output files.",
# )
# def featurize_seqs(
#     hk_df,
#     output
# ):
#     """ Primarily for PB/MB cells - add DMS features to the dataframe
#     including num mutations and additive Kd Scores for each sequence (row) in df
#
#     HK_DF: A dataframe with paired HK seqs
#     VARIANT_SCORES: Path to variant scores csv in DMS repository
#     NAIVE_SITES: Path to sites csv in DMS repository
#     """
#     # DMS single mutant scores
#     dms_df = pd.read_csv(
#         variant_scores, index_col="mutation", dtype=dict(position_IMGT=pd.Int16Dtype())
#     )
#     # remove linker sites
#     dms_df = dms_df[dms_df.chain != "link"]
#
#     # position maps for scFv
#     pos_df = pd.read_csv(
#         naive_sites,
#         dtype=dict(site=pd.Int16Dtype()),
#         index_col="site_scFv",
#     )
#     # position maps for heavy and light chain that can be used in gctree render
#     igh_pos_map = pos_df.loc[pos_df.chain
#                              == "H", "site"].reset_index(drop=True)
#     igk_pos_map = pos_df.loc[pos_df.chain
#                              == "L", "site"].reset_index(drop=True)
#
#     # TODO
#     pass


# WSD note: pending discussion, but I think we may not need this, so commenting out
##### MERGE FEATURIZED TREE SEQS
# @cli.command("merge-featurized-tree-seqs")
# @click.argument("hk_df", type=click.Path(exists=True))
# @click.argument("featurized_gctree_file", type=click.Path(exists=True))
# @click.option(
#     "--output",
#     "-o",
#     type=click.Path(exists=True),
#     default=".",
#     help="Path to write output files.",
# )
# def merge_featurized_tree_seqs(
#     hk_df,
#     featurized_gctree_file,
#     output
# ):
#     """This function should take the original HK df containing
#     extant seqs for any single gctree inference, and the rank1 featurized
#     tree created using `gctree infer` &  `featurize_nodes`,
#     It will then extract the DMS and Additive Kd info from the nodes,
#     merging all the information (including putative nodes) into
#     a final HK df containing all information
#
#     HK_DF: A dataframe with paired HK seqs
#     GCTREE_FILE: Path to pickled gctree.CollapsedTree object
#     """
#
#     # TODO
#     pass


# TODO We need this until we get setup.py
# for a real mf python package.
if __name__ == '__main__':
    cli()
