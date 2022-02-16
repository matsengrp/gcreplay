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
import click
import re

#################################
# HELPERS
#################################

# DEFine a function that collapses a numpy array, like the one below
def collapser()
    """
    it's got to take a 
    it's got to return a list
    lambda c:[c for c in partis_airr.columns if c not in ["counts", "rank", "num"]
             [np.max(x.Score),
                 df.loc[x.Score.idxmax(),'Element'],
                 df.loc[x.Score.idxmax(),'Case'],
                 np.min(x.Evaluation)])
    """
    
    # def apply(fn, x)

    #return (
    #    df.groupby('Group').apply(

    #      # lambda x: -> []

    #      .apply(pd.Series)
    #      .rename(columns={0:'Max_score_value',
    #                       1:'Max_score_element',
    #                   2:'Max_score_case',
    #                   3:'Min_evaluation'})
    #      .reset_index()
    #)


def parse_nextflow_header(header: str):
    """parse fasta header and return rank, counts, well, plate, and chain, 
    we do not expect the `>` to be included in the header"""
    if "unmatched" in header: return -1
    
    pr, date, plate, well, chain, _, rank_count = header.split(".")
    bcr_ranking, bcr_count  = rank_count.split("-")
    
    return {
        "identifier":header, 
        "plate":plate, 
        "well":well, 
        "chain":chain, 
        "rank":int(bcr_ranking), 
        "count":int(bcr_count)
    }


# probably dont need this? 
def parse_tatsuya_header(header: str):
    """parse fasta header and return rank, counts, well, plate, and chain, 
    we do not expect the `>` to be included in the header"""
    if "unmatched" in header: return -1
    
    prefix, rank_count = header.split("_")
    bcr_ranking, bcr_count  = rank_count.split("-")
    
    # pattern match for plate, well, chain
    plate = re.search(r'([P][0-9]+)', prefix).group(0)
    well = re.search(r'([A-L][0-9]+)', prefix).group(0)
    chain = prefix[-1]
    return {
        "identifier":header, 
        "plate":plate, 
        "well":well, 
        "chain":chain, 
        "rank":int(bcr_ranking), 
        "count":int(bcr_count)
    }


def bcr_fasta_to_df(fasta_fp, id_parse_fn, **kwargs)    
    """convert a fasta file pointer to dataframe after parsing the id with some function returning
    the columns defined (less the sequence column)"""
    
    columns = ['identifier', 'plate', 'well', 'chain', 'rank', 'count', 'sequence']
    ret  = pd.DataFrame({c:[] for c in columns})
    with open(fasta_fp) as fasta_file:  # Will close handle cleanly
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            bcr_meta = id_parse_fn(seq_record.id)
            if bcr_meta != -1:
                bcr_meta["sequence"] = str(seq_record.seq)
                ret = ret.append(pd.Series(bcr_meta), ignore_index=True)
    return ret


def plot_venn_stub(
    df1:pd.DataFrame, 
    df2:pd.DataFrame, 
    feature_groups:list, 
    out="venn.png"
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
    This function should take in a collapsed fasta, parse it into a dataframe,
    and output a fasta 
    """
    pass




#################################
# CLI
#################################

"""
$ python hello.py --count=3
Your name: John
Hello John!
Hello John!
Hello John!
"""

@click.command()
@click.option('--count', default=1, help='Number of greetings.')
@click.option('--name', prompt='Your name',
              help='The person to greet.')
def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for x in range(count):
        click.echo(f"Hello {name}!")

if __name__ == '__main__':
    hello()





