#!/usr/bin/env python
"""
@file: gcreplay-tools.py
"""

import os
import glob

import pandas as pd
import click

##### MERGE FEATURIZED TREE SEQS
@click.command()
@click.option("--glob-pattern", type=str, default="*/")
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
def merge_results(glob_pattern, output_gcdf, output_gctree):
    """
    merge all the gc results.
    """

    gcdf, gctree = pd.DataFrame(), pd.DataFrame()
    for gct_out in glob.glob(glob_pattern):
        gcdf = pd.concat([gcdf, pd.read_csv(f"{gct_out}/observed_seqs.csv")])

        # loop through each of the ranking coeff sub-directories under the gct_out directory
        uid = gct_out.split("/")[0]
        for rc in glob.glob(f"{gct_out}/*"):
            if os.path.isdir(rc):
                ranking_coeff_strategy = rc.split("/")[-1]
                node_data = pd.read_csv(f"{rc}/node_data.csv")
                node_data.insert(0, "uid", uid)
                node_data.insert(1, "ranking_coeff_strategy", ranking_coeff_strategy)
                gctree = pd.concat([gctree, node_data])

    gcdf.to_csv(output_gcdf, index=False)
    gctree.to_csv(output_gctree, index=False)

if __name__ == "__main__":
    merge_results()
