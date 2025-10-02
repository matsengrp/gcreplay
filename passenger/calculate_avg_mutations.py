"""
Calculate simple average mutations per sequence for IGH and IGK experiments.

This script computes the total number of mutations observed divided by the
number of sequences for each experiment, providing a straightforward measure
of mutation load.
"""

import pandas as pd
import numpy as np

print("=" * 70)
print("IGH: Heavy Chain Analysis")
print("=" * 70)

# Load IGH data
counts_by_base_igh = pd.read_csv('output/igh_counts_by_base.csv')
read_counts_igh = pd.read_csv('output/igh_read_counts.csv')

# Group by nickname and sum
by_base_dfs_igh = dict(list(counts_by_base_igh.groupby("nickname")))
by_site_igh = {}
for name, df in by_base_dfs_igh.items():
    df.drop(columns=["nickname"], inplace=True)
    df.set_index("position", inplace=True)
    df.replace(-1, float("nan"), inplace=True)
    # Sum across bases (A, C, G, T) for each position
    by_site_igh[name] = df.sum(axis=1, skipna=True)

by_site_df_igh = pd.DataFrame(by_site_igh)

# Calculate totals before any pseudocount
total_mutations_igh = by_site_df_igh.sum(axis=0)
total_sequences_igh = read_counts_igh['read_count'].values

results_igh = pd.DataFrame({
    'nickname': read_counts_igh['nickname'],
    'total_mutations': total_mutations_igh.values,
    'num_sequences': total_sequences_igh,
    'avg_mutations_per_seq': total_mutations_igh.values / total_sequences_igh
})

print("\nPer-experiment results:")
print(results_igh.to_string(index=False))

print(f"\n{'=' * 70}")
print(f"Overall IGH Summary:")
print(f"{'=' * 70}")
print(f"Total mutations (all experiments): {total_mutations_igh.sum():.0f}")
print(f"Total sequences (all experiments): {total_sequences_igh.sum()}")
print(f"Overall avg mutations per sequence: {total_mutations_igh.sum() / total_sequences_igh.sum():.4f}")

print("\n\n")
print("=" * 70)
print("IGK: Kappa Light Chain Analysis")
print("=" * 70)

# Load IGK data
prefixes = ["8a", "9a", "11a"]
read_counts_igk = np.array([8662, 45079, 74396])  # From the notebook

by_site_igk = {}
for prefix in prefixes:
    df = pd.read_csv(f'output/{prefix}_counts_by_base.csv', index_col=0)
    df = df.replace(-1, float("nan"))
    by_site_igk[prefix] = df.sum(axis=1, skipna=True)

by_site_df_igk = pd.DataFrame(by_site_igk)

# Calculate totals before any pseudocount
total_mutations_igk = by_site_df_igk.sum(axis=0)

results_igk = pd.DataFrame({
    'nickname': prefixes,
    'total_mutations': total_mutations_igk.values,
    'num_sequences': read_counts_igk,
    'avg_mutations_per_seq': total_mutations_igk.values / read_counts_igk
})

print("\nPer-experiment results:")
print(results_igk.to_string(index=False))

print(f"\n{'=' * 70}")
print(f"Overall IGK Summary:")
print(f"{'=' * 70}")
print(f"Total mutations (all experiments): {total_mutations_igk.sum():.0f}")
print(f"Total sequences (all experiments): {read_counts_igk.sum()}")
print(f"Overall avg mutations per sequence: {total_mutations_igk.sum() / read_counts_igk.sum():.4f}")
