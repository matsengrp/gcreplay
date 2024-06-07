"""Parse metadata files."""

import pandas as pd
import pathlib

path = pathlib.Path(__file__).parent
df = pd.concat([pd.read_csv(path / f"../../metadata.PR{pr}.csv", dtype=str)
                for pr in (1, 2)], ignore_index=True)

df.rename(columns={"imm_duration": "time", "gc": "GC"}, inplace=True)

df.sort_values(["cell_type", "time", "PR", "mouse", "GC"], inplace=True)

df.index = pd.Index("PR" + df.PR + "_mouse" + df.mouse + "_GC" + df.GC, name="uid")
assert df.index.is_unique

# renamed mice and GCs
df_renamed = df.copy()
# Ashni says "The mouse numbers in 2.01+2.02 refer to the same mice and those in 2.03+2.04 refer to the same mice"
def pr_group(pr):
    if pr.startswith("1."):
        return "1"
    if pr in ("2.01", "2.02"):
        return "2.01/2.02"
    if pr in ("2.03", "2.04"):
        return "2.03/2.04"
    return pr
_pr_groups = df_renamed.PR.apply(pr_group)
_mice = _pr_groups + "_" + df_renamed.mouse
_gc = df_renamed.PR + "_" + df_renamed.mouse + "_" + df_renamed.GC
df_renamed.mouse = (pd.factorize(_mice)[0] + 1).astype(str)
df_renamed.GC = (pd.factorize(_gc)[0] + 1).astype(str)
assert not df_renamed[['mouse', 'GC']].duplicated().any()


# df.query("(strain == 'wt') & (cell_type == 'GC') & (time != 'w10')", inplace=True)


def main():
    df_renamed.to_csv(sys.stdout)


if __name__ == "__main__":
    import sys
    main()
