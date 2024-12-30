#!/usr/bin/env python
"""
@file: wrangle-annotation.py
"""

import pandas as pd
import click
from click import Path

from utils import partis_airr_to_drop, final_HK_col_order
from utils import infer_igh_isotypes, bcr_fasta_to_df, parse_nextflow_header


@click.command()
@click.option(
    "--igh-airr",
    type=Path(exists=True),
    required=True,
    help="igh airr output from partis",
)
@click.option(
    "--igk-airr",
    type=Path(exists=True),
    required=True,
    help="igk airr output from partis",
)
@click.option(
    "--input-fasta",
    type=Path(exists=True),
    required=True,
    help="the fasta fed into partis for annotation",
)
@click.option(
    "--gc-metadata",
    type=Path(exists=True),
    required=True,
    help="the key file for merging heavy and light chains",
)
@click.option(
    "--ngs-id",
    required=True,
    help="the ngs id to wrangle",
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Path to write the fasta.",
)
def wrangle_annotation(igh_airr, igk_airr, input_fasta, gc_metadata, ngs_id, output):
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
        f"{r.ngs_date}{r.plate}{r.well}{r.chain}" for i, r in partis_airr.iterrows()
    ]

    partis_airr = partis_airr.rename(
        {
            "v_call": "V",
            "d_call": "D",
            "j_call": "J",
            "productive": "Productive",
            "junction_aa": "AAjunction",
            "sequence_id": "fasta_header",
            "plate": "miseq_plate",
            "sequence": "partis_sequence",
        },
        axis=1,
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
        infer_igh_isotypes(row.fasta_seq, 3) if row.locus == "IGH" else "IgK"
        for idx, row in partis_airr.iterrows()
    ]
    no_iso = partis_airr.loc[partis_airr.isotype.isnull()].index
    partis_airr.drop(no_iso, inplace=True)

    # Some cells are spread across wells.
    # This means we need to merge the counts across plates
    # where wells are identical
    assert len(partis_airr) == len(set(partis_airr.index))

    gc_metadata_df = pd.read_csv(gc_metadata).query("ngs_id == @ngs_id")
    for idx, row in gc_metadata_df.iterrows():
        for chain, column in zip(["IGH", "IGK"], ["hc_barcode", "lc_barcode"]):
            key_chain_barcodes = [int(c) for c in str(row[column]).split(".")]

            if len(key_chain_barcodes) > 1:

                # Let's make sure this code isn't effecting any unwanted
                # normal GC's. I assume it won't
                assert row.gc == 53

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
                    sorted_sum_well_seqs = well_summed_seqs.sort_values(
                        "counts", ascending=False
                    )
                    winning_seq = sorted_sum_well_seqs.index.values[0]
                    winning_counts = sorted_sum_well_seqs.counts.values[0]
                    winning_original_entries = well_multiplate_df.query(
                        f"fasta_seq == '{winning_seq}'"
                    )
                    winning_index = winning_original_entries.index.values[0]
                    to_throw.remove(winning_index)
                    partis_airr.loc[winning_index, "counts"] = winning_counts
                    partis_airr.loc[winning_index, "rank"] = 1

                # now drop all the non winning BCR's
                partis_airr.drop(to_throw, axis=0, inplace=True)

    partis_airr.query("(rank == 1)", engine="python", inplace=True)
    for well, well_df in partis_airr.groupby("ID"):
        assert len(well_df) == 1

    # Merge heavy and light chain

    GC_df = pd.DataFrame()
    for idx, row in gc_metadata_df.iterrows():
        key_row = row.row.split(".")
        key_col = [int(c) for c in row.col.split(".")]
        key_hc_barcodes = [int(c) for c in str(row.hc_barcode).split(".")]
        key_lc_barcodes = [int(c) for c in str(row.lc_barcode).split(".")]

        HC = partis_airr.query(
            (
                f"(locus == 'IGH') & "
                f"(barcode.isin({key_hc_barcodes})) & "
                f"(column.isin({key_col})) & "
                f"(row.isin({key_row}))"
            ),
            engine="python"
        )
        LC = partis_airr.query(
            (
                f"(locus == 'IGK') & "
                f"(barcode.isin({key_lc_barcodes})) & "
                f"(column.isin({key_col})) & "
                f"(row.isin({key_row}))"
            ),
            engine="python"
        )

        GC = HC.merge(LC, on="well", suffixes=("_HC", "_LC"))
        GC["uid"] = row.uid
        GC["ngs_id"] = row.ngs_id
        GC["mouse"] = row.mouse
        GC["gc"] = row.gc
        GC["node"] = row.node
        GC["cell_type"] = row.cell_type
        GC["imm_duration"] = row.imm_duration

        GC_df = pd.concat([GC_df, GC])

    GC_df.loc[:, "ID_HK"] = [f"{i}K" for i in GC_df["ID_HC"]]
    GC_df = GC_df.query("chain_HC != chain_LC")

    # TODO This should be somewhere else, but it's a qc on number of Ns in a sequence.
    throw_seq = [
        idx
        for idx, row in GC_df.iterrows()
        if (row.seq_nt_HC + row.seq_nt_LC).count("n") > 30
    ]
    GC_df.drop(throw_seq, axis=0, inplace=True)

    ret = GC_df.loc[:, [c for c in final_HK_col_order if c in GC_df.columns]]
    ret.to_csv(output, index=False)

if __name__ == "__main__":
    wrangle_annotation()
