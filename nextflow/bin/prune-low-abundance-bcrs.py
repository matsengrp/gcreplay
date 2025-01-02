#!/usr/bin/env python
"""
@file: prune_low_abundance_bcrs.py
"""


from Bio import SeqIO
import click
from click import Path


@click.command()
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
def prune_low_abundance_bcrs(fasta: str, count_threshold: int, output: str):
    """Keep only the sequences with counts above the threshold."""

    with open(fasta, "r") as f1, open(output, "w") as f2:
        for seq_record in SeqIO.parse(f1, 'fasta'):  # (generator)
            rank, count = str(seq_record.id).split("-")
            if int(count) < count_threshold: break
            f2.write(f">{seq_record.id}\n{seq_record.seq}\n")

if __name__ == '__main__':
    prune_low_abundance_bcrs()