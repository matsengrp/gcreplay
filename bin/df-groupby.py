#!/usr/bin/env python
"""
@file: df-groupby.py
"""

import pandas as pd
import click
from click import Path

@click.command()
@click.option(
    "--dataframe",
    "-df",
    type=Path(exists=True),
    required=True,
    help="dataframe (csv) to query",
)
@click.option(
    "--column",
    "-c",
    type=str,
    required=True,
    default="uid",
    help="the key file for merging heavy and light chains",
)
@click.option("--sample", type=int, default=None)
def df_groupby(dataframe, column, sample):
    """
    Simply, a CLI wrapper for pandas DataFrame.groupby on certain columns.
    """
    df = pd.read_csv(dataframe)
    for group, groupdf in df.groupby(column):
        if sample is not None:
            groupdf.sort_values(by="ID_HK").sample(
                min(10, len(groupdf)), random_state=1
            ).to_csv(f"{group}.csv", index=False)
        else:
            groupdf.sort_values(by="ID_HK").to_csv(
                f"{group}.csv", index=False
            )

if __name__ == "__main__":
    df_groupby()