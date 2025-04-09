import argparse
import pandas as pd
import os
from utils import naive_hk_bcr_nt



parser = argparse.ArgumentParser(
    description='A script that takes in an observed-seqs csv and creates fasta for input to beast.nf'
)

def gc_df_to_fasta(
        gc_hk_df, 
        output, 
        header_col=["ID_HK"], 
        sequence_col=["seq_nt_HC", "seq_nt_LC"],
        add_naive=True
):
    """
    A function to concatinate specified columns from
    the germinal center dataframe (output by
    `wrangle_partis_annotation` command) for both headers
    and sequences.
    """

    # gather the concatinated columns for headers
    headers = gc_hk_df[header_col].apply(lambda x: "".join(x), axis=1)
    sequences = gc_hk_df[sequence_col].apply(
        lambda x: "".join(x), axis=1)

    # write fasta
    with open(output, "w") as fasta:
        if add_naive:
            fasta.write(f">naive\n{naive_hk_bcr_nt}\n")
        for header, sequence in zip(headers, sequences):
            fasta.write(f">{header}\n{sequence}\n")


# Adding string arguments
parser.add_argument('--df', type=str, help='')
parser.add_argument('--output_dir', type=str, help='')

if __name__ == "__main__":

    args = parser.parse_args()
    df = pd.read_csv(args.df)
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    for uid, sample_df in df.groupby("uid"):
        output = os.path.join(args.output_dir, f"{uid}.fasta")
        gc_df_to_fasta(
            sample_df,
            output
        )


