
#!/usr/bin/env python
"""
@file: utils.py

"""

import warnings
import os
import sys
import re
import pickle
import glob
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
# import click
# from click import Choice, Path, command, group, option, argument
import regex

import matplotlib
#matplotlib.use("Qt5Agg")


##################
# CONSTANTS
##################

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
    "aa_substitutions_IMGT",

    "delta_bind",
    "delta_expr"

    "n_nt_mutations_HC", "n_nt_mutations_LC", "IgH_nt_mutations", "IgK_nt_mutations",
    "n_aa_mutations_HC", "n_aa_mutations_LC", "IgH_aa_mutations", "IgK_aa_mutations",

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
]

#########
# HELPERS
#########

def fasta_to_df(f):
    """simply convert a fasta to dataframe"""

    ids, seqs = [], []
    with open(f) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            ids.append(seq_record.id)
            seqs.append(str(seq_record.seq))
    return pd.DataFrame({"id":ids, "seq":seqs})


def bcr_fasta_to_df(fasta_fp, id_parse_fn, **kwargs):
    """convert a fasta file pointer to dataframe after
    parsing the id with some function returning
    the columns defined (less the sequence column)"""

    columns = [
        'sequence_id',
        'ngs_date',
        'plate',
        'barcode',
        'well',
        'row',
        'column',
        'chain',
        'rank',
        'counts',
        'fasta_seq'
    ]

    ret = {c: [] for c in columns}
    with open(fasta_fp) as fasta_file:  # Will close handle cleanly
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            bcr_meta = id_parse_fn(seq_record.id)
            if bcr_meta != -1:
                bcr_meta["fasta_seq"] = str(seq_record.seq)
                for key, item in bcr_meta.items():
                    ret[key].append(item)

    return pd.DataFrame(ret)


def parse_nextflow_header(header: str):
    """parse fasta header and return rank, counts, well, plate, and chain,
    we do not expect the `>` to be included in the header"""
    if "unmatched" in header:
        return -1

    pr, date, plate, well, chain, _, rank_count = header.split(".")
    bcr_ranking, bcr_count = rank_count.split("-")

    return {
        "sequence_id": header,
        "ngs_date": date,
        "plate": plate,
        "barcode": int(plate[1:]),
        "well": well,
        "row": well[0],
        "column": int(well[1:]),
        "chain": chain,
        "rank": int(bcr_ranking),
        "counts": int(bcr_count)
    }



def infer_igh_isotypes(bcr_sequence, num_mm=2):
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

    assert len(naive_aa) == len(aa)
    return [
        f"{aa1}{pos_map[pos]}{chain_annotation}{aa2}"
        for pos, (aa1, aa2) in enumerate(zip(naive_aa, aa))
        if aa1 != aa2
    ]

def nt_mutations(naive_nt, nt, site_idx_offset=0):
    """Nucleotide substitutions between two sequences, in sequential coordinates."""
    assert len(naive_nt) == len(nt)
    return [
        f"{nt1}{pos + site_idx_offset}{nt2}"
        for pos, (nt1, nt2) in enumerate(zip(naive_nt, nt))
        if nt1 != nt2
    ]


# move this to annotation as well
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
        key_hc_barcodes = [int(c) for c in str(row.hc_barcode).split(".")]
        key_lc_barcodes = [int(c) for c in str(row.lc_barcode).split(".")]

        HC = cell_df.query(
            (
                f"(locus == 'IGH') & "
                f"(barcode.isin({key_hc_barcodes})) & "
                f"(column.isin({key_col})) & "
                f"(row.isin({key_row}))"
            ),
            engine="python"
        )
        LC = cell_df.query(
            (
                f"(locus == 'IGK') & "
                f"(barcode.isin({key_lc_barcodes})) & "
                f"(column.isin({key_col})) & "
                f"(row.isin({key_row}))"
            ),
            engine="python"
        )

        GC = HC.merge(LC, on="well", suffixes=("_HC", "_LC"))
        GC["HK_key_mouse"] = row.mouse
        GC["HK_key_gc"] = row.gc
        GC["HK_key_node"] = row.node
        GC["HK_key_cell_type"] = row.cell_type
        GC["HK_key_plate"] = row.plate

        GC_df = pd.concat([GC_df, GC])

    GC_df.loc[:, "ID_HK"] = [f"{i}K" for i in GC_df["ID_HC"]]
    GC_df = GC_df.query("chain_HC != chain_LC")

    return GC_df



