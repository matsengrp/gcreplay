#!/usr/bin/env python
"""
@file: gcreplay-tools.py
"""

from utils import fasta_to_df
import click
from click import Path


@click.command()
@click.option(
    "--fasta",
    type=Path(exists=True),
    required=True,
    help="a sequence count rank collapsed fasta with header format <rank>-<count>",
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Path to write the fasta.",
)
def sort_fasta(fasta: str, output: str):
    """sort a fasta file by sorting the headers using
    pandas lexiconographical sort function"""

    fasta_df = fasta_to_df(fasta).sort_values(by="id")
    with open(output, "w") as f2:
        for idx, row in fasta_df.iterrows():
            f2.write(f">{row.id}\n{row.seq}\n")

if __name__ == '__main__':
    sort_fasta()