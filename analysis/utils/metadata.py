"""Parse metadata files."""

import numpy as np
import pandas as pd


df = pd.concat([pd.read_csv(f"../metadata.PR{pr}.csv", dtype=str)
                      for pr in (1, 2)], ignore_index=True)

df.rename(columns={"imm_duration": "time", "gc": "GC"}, inplace=True)

df.sort_values(["cell_type", "time", "PR", "mouse", "GC"], inplace=True)

df.index = "PR" + df.PR + "_mouse" + df.mouse + "_GC" + df.GC
assert df.index.is_unique

# renamed mice and GCs
df_renamed = df.copy()
_mice = df_renamed.PR + "_" + df_renamed.mouse
_gc = df_renamed.PR + "_" + df_renamed.mouse + "_" + df_renamed.GC
df_renamed.mouse = (pd.factorize(_mice)[0] + 1).astype(str)
df_renamed.GC = (pd.factorize(_gc)[0] + 1).astype(str)


# df.query("(strain == 'wt') & (cell_type == 'GC') & (time != 'w10')", inplace=True)
