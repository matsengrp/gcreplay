"""
bug tracking
"""

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import glob
import os
import sys
#import gcreplay_tools as gct

run1 = sys.argv[1]
run2 = sys.argv[2]

print("\n\n############################################")
print("POST WRANGLE DIFFS")
total = pd.DataFrame()
for f in glob.glob(f"{run1}/single_gc_wrangle/*"):

    fn = os.path.basename(f)
    print(f)
    r1 = pd.read_csv(f"{run1}/single_gc_wrangle/{fn}", index_col = "ID_HK")
    r2 = pd.read_csv(f"{run2}/single_gc_wrangle/{fn}", index_col = "ID_HK")
    diff = pd.concat([r1, r2]).drop_duplicates(keep=False)
    if len(diff) == 0: continue
    total = pd.concat([total,diff])
    print(f"there a {len(diff)} different rows between the two runs, they are\n{diff}")
    print("=================================")

print(total)
total.to_csv("total_PR1.6_diffs_HK_post_wrangle.csv")

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

print("\n\n############################################")
print("PARTIS INPUT FASTA DIFFS")

def fasta_to_df(f):
    ids, seqs = [], []
    with open(f) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            ids.append(seq_record.id)
            seqs.append(str(seq_record.seq))
    return pd.DataFrame({"id":ids, "seq":seqs})

r1_input_fasta = fasta_to_df(f"{run1}/partis_annotation/PR-1-6.fasta") 
r2_input_fasta = fasta_to_df(f"{run2}/partis_annotation/PR-1-6.fasta") 
assert len(r1_input_fasta) == len(set(r1_input_fasta["id"]))

diff = pd.concat([r1_input_fasta, r2_input_fasta]).drop_duplicates(keep=False)
if len(diff) == 0: print("No differences in input fasta")


print("\n\n############################################\n\n")
print("PARTIS OUTPUT DIFFS (single chain tsv outputs)")

for chain in ["igk", "igh"]:

    print(chain, "chain output differences")
    print("**********************************************")
    r1 = pd.read_csv(f"{run1}/partis_annotation/PR-1-6/engrd/single-chain/partition-{chain}.tsv", index_col = "sequence_id", sep="\t")
    r2 = pd.read_csv(f"{run2}/partis_annotation/PR-1-6/engrd/single-chain/partition-{chain}.tsv", index_col = "sequence_id", sep="\t")
    diff = pd.concat([r1, r2]).drop_duplicates(keep=False)
    if len(diff) == 0: continue
    diff.to_csv(f"total_PR1.6_{chain}_diffs_partis_anno.csv")
    
    for g, gdf in diff.groupby("sequence_id"):
        pgdf = gdf[["sequence", "seqs_aa", "n_mutations"]]
        input_seq = r1_input_fasta[r1_input_fasta["id"]==g]["seq"].values[0]
        print("DIFFERENCE")
        print(f"input fasta id: {g}")
        print(f"input fasta seq: {input_seq}")
        print(f"difference between annotations:\n{pgdf}") 
        print("=================================")

    


    

#r1_wrangle = gct.wrangle_annotation(
#    igh_airr = "wrangle-1/partis_annotation/PR-1-6/engrd/single-chain/partition-igh.tsv",
#    igk_airr = "wrangle-1/partis_annotation/PR-1-6/engrd/single-chain/partition-igk.tsv",
#    input_fasta = "wrangle-1/partis_annotation/PR-1-6.fasta",
#    key_file = "PR1_6-key.csv",
#    output = "test.csv"
#)
#
#r2_wrangle = gct.wrangle_annotation(
#    igh_airr = "wrangle-2/partis_annotation/PR-1-6/engrd/single-chain/partition-igh.tsv",
#    igk_airr = "wrangle-2/partis_annotation/PR-1-6/engrd/single-chain/partition-igk.tsv",
#    input_fasta = "wrangle-2/partis_annotation/PR-1-6.fasta", 
#    key_file = "PR1_6-key.csv",
#    output = "test.csv"
#)

#diff = pd.concat([r1_wrangle, r2_wrangle]).drop_duplicates(keep=False)
#print(f"there a {len(diff)} different rows between the two wrangles, they are\n{diff}")
#diff.to_csv()
##print(diff.seq_aa.unique())

