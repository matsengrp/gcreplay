"""
@file: gcreplay-tools

Right now a sandbox for throwing code into
it'll be a simple click CLI with all the logic right here. 
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import *
from matplotlib_venn import venn2, venn3_circles
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
import click
from click import Choice, Path, command, group, option, argument
import regex
import re


#################################
# GLOBALS
#################################

# TODO we should probably just store the important sequences in fasta
naive_hk_bcr_nt = "GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCT \
CAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGA \
TCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTA \
CTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTAC \
CTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATG \
TCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAGACATTGTGATGACtCAGTCTCAAAAATT \
CATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACT \
AATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCT \
ACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCAC \
CATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCT \
CTCACGTTCGGCTCGGGGACtAAGCTaGAAATAAAA"

naive_hk_bcr_aa = "EVQLQESGPSLVKPSQTLSLTCSVTGDSITSGYWNWIRKFPGNKLEYMGYISYSG \
STYYNPSLKSRISITRDTSKNQYYLQLNSVTTEDTATYYCARDFDVWGAGTTVTVSSDIVMTQS \
QKFMSTSVGDRVSVTCKASQNVGTNVAWYQQKPGQSPKALIYSASYRYSGVPDRFTGSGSGTDF \
TLTISNVQSEDLAEYFCQQYNSYPLTFGSGTKLEIK"

partis_airr_to_drop = [
    "j_sequence_end",
    "j_germline_start",
    "j_germline_end",  
    "rev_comp",
    "junction",
    "clone_id",
    "vj_in_frame",
    "stop_codon",
    "np1",
    "np2",
    "duplicate_count",
    "cdr3_start",
    "cdr3_end",
    "cell_id",
    "v_support",
    "v_identity",
    "v_sequence_start",
    "v_sequence_end",
    "v_germline_start",
    "v_germline_end",
    "d_support",
    "d_identity",
    "d_sequence_start",
    "d_sequence_end",
    "d_germline_start",
    "d_germline_end",
    "j_support",
    "j_identity",
    "j_sequence_start",
    "j_sequence_end",
    "j_germline_start",
    "j_germline_end",
]

# TODO
isotype_motifs = {
    
}

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

    ret  = pd.DataFrame({c:[] for c in columns})
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
    if "unmatched" in header: return -1
    
    pr, date, plate, well, chain, _, rank_count = header.split(".")
    bcr_ranking, bcr_count  = rank_count.split("-")
    
    return {
        "sequence_id":header, 
        "date":date,
        "plate":plate, 
        "barcode":int(plate[1:]),
        "well":well, 
        "row": well[0],
        "column": int(well[1:]),
        "chain":chain, 
        "rank":int(bcr_ranking), 
        "counts":int(bcr_count)
    }


def plot_venn_stub(
    df1:pd.DataFrame, 
    df2:pd.DataFrame, 
    feature_groups:list, 
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


def threshold_fasta_sequence_abundance(fasta):
    """
    str -> str

    This function should take in a collapsed fasta, parse it into a dataframe,
    and output a fasta 
    """
    pass


# TODO
def infer_igh_isotypes(bcr_sequence, num_mm=1):
    """Use fuzzy string matching to infer the isotype of bcr sequences
    in the cell database - allowing for `num_mm` mismatches. This function
    will modify the passed in dataframe and will return the first match it finds
    """
    
    motif_map = {
        "IgM": ["atgtcttccccct"],
        "IgG1" : ["atggtgaccctggg"],
        "IgG2" : ["ggctcctcggtgactcta", "tcggtgactctagg", "gtgtggagatacaactgg"],
        "igG3" : ["cccttggtccctggctgcggtgacacat", "cacatctggatcctcggtgaca"],
        "IgA" : ["tctgcgagaaatcccaccatcta"]
    }

    # match and return first valid isotype
    for key, value in motif_map.items():
        for motif in value:
            m = regex.findall("({m}){{e<={nm}}}".format(m=motif, nm=num_mm), bcr_sequence.lower())
            if len(m) > 0:
                return key
    # TODO Log
    print(f"WARNING: no isotype match for: {bcr_sequence}")

    # if no istype motif is found
    return np.nan


def _test_infer_igh_isotypes() -> None:
    assert infer_igh_isotypes("cctgaagtctgcgagaaatcccaccatctatctta") == "IgA"
    assert infer_igh_isotypes("cctgaagtctgcgagatatcccaccatctatctta") == "IgA"


# TODO
def infer_mutations(cell_data: pd.DataFrame) -> None:
    """
    """
    pass


# TODO
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
        #print(f"(locus == 'IGH') & (barcode == {row.hc_barcode}) & (column.isin({key_col})) & (row.isin({key_row}))")

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

        GC_df = GC_df.append(GC)
   
    # TODO we're missing just to compared to tatsuya, I believe - check this

    GC_df.loc[:, "ID_HK"] = [f"{i}K" for i in GC_df["ID_HC"]]

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
    a gcreplay experiment.
    """
    pass

# TODO option, filter unproductive cells
# TODO option, mm in isotype inference
# TODO option, 
@cli.command()
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
def wrangle_annotation(
    igh_airr, igk_airr, input_fasta, key_file
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
    partis_airr.loc[:, "seq_nt"] = [seq.lower() for seq in partis_airr["sequence"]]

    # drop unnecessary columns
    partis_airr.drop(partis_airr_to_drop, axis=1, inplace=True)


    # parse the fasta to merge in the important information 
    # TODO add plate int to this as 'barcode', keep plate for the hell of it
    parsed_input_fasta = bcr_fasta_to_df(input_fasta, parse_nextflow_header)
    

    # Merge in the parsed fasta columns 
    # TODO once we're getting rid of unmatched we can do a left join? assert this.
    partis_airr = partis_airr.merge(parsed_input_fasta, how="inner", on="sequence_id")
    partis_airr["ID"] = [
        f"{r.date}{r.plate}{r.well}{r.chain}"
        for i, r in partis_airr.iterrows()
    ]


    partis_airr = partis_airr.rename(
        {
            "v_call":"V",
            "d_call":"D",
            "j_call":"J",
            "productive":"Productive",
            "junction_aa":"AAjunction",
            "sequence_id" : "fasta_header"
        },
        axis=1
    )

    # compute lengths
    partis_airr.loc[:, "seq_nt_length"] = [len(seq) for seq in partis_airr["seq_nt"]]
    partis_airr.loc[:, "seq_aa_length"] = [len(seq) for seq in partis_airr["seq_aa"]]

    query_string = "(locus == 'IGH' & seq_nt_length == 337) | (locus == 'IGK' & seq_nt_length == 321)"

    # Trim the rest of partis IGH
    partis_airr.query(query_string, engine="python", inplace=True)
    seq_nts = []
    for idx, row in partis_airr.iterrows():
        seq_nts.append(row.seq_nt[:-1] if row.locus == "IGH" else row.seq_nt)
    partis_airr.loc[:, "seq_nt"] = seq_nts

    # re-compute lengths
    partis_airr.loc[:, "seq_nt_length"] = [len(seq) for seq in partis_airr["seq_nt"]]
    partis_airr.loc[:, "seq_aa_length"] = [len(seq) for seq in partis_airr["seq_aa"]]

    partis_airr.loc[:, "isotype"] = [
        infer_igh_isotypes(row.seq_input) 
        if row.locus == "IGH" else "IgK"
        for idx, row in partis_airr.iterrows()
    ]
   
    # TODO shall we export those which are not rank 1, count 10?
    partis_airr.query("(rank == 1) & (counts >= 10)", engine="python", inplace=True)
    for well, well_df in partis_airr.groupby("ID"):
        assert len(well_df) == 1

    GC_df = merge_heavy_light_chains(partis_airr, pd.read_csv(key_file)
    

