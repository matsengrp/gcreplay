# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # Spatial correlation of passenger-allele mutations
#
# First-cut analysis for
# [matsengrp/gcreplay#179](https://github.com/matsengrp/gcreplay/issues/179):
# compute and plot the spatial correlation function
# $f(r) = n(r)/n_m(r)$ on gcreplay passenger (chigy stop allele) reads,
# where $n(r)$ is the observed number of mutation pairs at genomic
# distance $r$ and $n_m(r)$ is the expected number under the
# independent-sites null $n_m(r) = N_\text{reads} \sum_{|i-j|=r,\,i<j} p_i p_j$.
#
# Because the passenger allele is out-of-frame and carries a stop,
# every read is selection-free, so any short-range enrichment of
# $f(r)$ above 1 cannot be attributed to epistasis and is a direct
# probe of the SHM machinery. This follows
# [Spisak et al. 2020](https://doi.org/10.1093/nar/gkaa825), which
# found up to 5× enrichment at short $r$ in productive and
# non-productive repertoire data, and complements the "multihit
# correction" in [Matsen et al. 2025 (DASM)](https://doi.org/10.1101/2025.10.21.683652)
# which lumps this phenomenon into three scalar rate multipliers
# per codon.
#
# The final target of this notebook is an HTML render (not committed);
# the notebook itself is also not committed.

# %% [markdown]
# ## Setup

# %%
import sys
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path.cwd()))
import passenger as P
from passenger import (
    CHIGY_HC,
    CHIGY_LC,
    blast_df_of_blast_files,
    mutation_frequency_by_position_of,
)

RNG = np.random.default_rng(20250409)
sns.set_theme(context="notebook", style="whitegrid")


# %% [markdown]
# ## Shared filter matching the per-dataset notebooks
#
# Every `passenger-{igh,8a,9a,11a}.ipynb` notebook applies the same
# `chigy_believable` filter:
#
# - `mutation_count < 10`
# - `n_count < 10`
# - `bookended_by_gaps == True`
# - `gap_segment_count ∈ {2, 3}`
#
# The `mutation_count < 10` cap is a property of the upstream pipeline;
# the issue's proposed "≤15 HC / ≤10 LC" all-load stratum is
# unreachable as a result, so we use `≤ 9` (all-load) and `≤ 5`
# (low-load, closer to single-SHM-event regime) strata.

# %%
MAX_MUTATION_COUNT = 10  # matches every per-dataset notebook
MAX_N_COUNT = 10
ALLOWABLE_GAP_SEGMENT_COUNTS = [2, 3]
LOW_LOAD_CAP = 5  # mutation_count <= 5


def apply_believable_filter(processed_stop_df: pd.DataFrame) -> pd.DataFrame:
    df = processed_stop_df[
        (processed_stop_df["mutation_count"] < MAX_MUTATION_COUNT)
        & (processed_stop_df["n_count"] < MAX_N_COUNT)
        & (processed_stop_df["bookended_by_gaps"] == True)
    ].copy()
    df = df[df["gap_segment_count"].isin(ALLOWABLE_GAP_SEGMENT_COUNTS)]
    return df


# %% [markdown]
# ## Load and align reads
#
# HC: every `passenger-blast/outname_[5-7]*.tsv`, aligned to
# `CHIGY_HC`, `query_length = 20` (matches `passenger-igh.ipynb`).
#
# LC: the three `8A`, `9A`, `11A` files, aligned to `CHIGY_LC` at the
# query lengths used in each per-dataset notebook (8A: 32, 9A: 30,
# 11A: 32).
#
# This re-runs `Passenger.processed_stop_df_of_blast_df` on the BLAST
# outputs already on disk, because per-read `mutation_positions` are
# not serialized by the existing aggregate notebooks.

# %%
HC_BLAST_PATHS = sorted(glob("passenger-blast/outname_[5-7]*.tsv"))
HC_BLAST_PATHS

# %%
LC_SPECS = [
    {"path": "passenger-blast/outname_AV_IgYKSTOP_8A_S2_R1_001_atleast-2.blast.tsv", "query_length": 32},
    {"path": "passenger-blast/outname_AV_IgYKSTOP_9A_S3_R1_001_atleast-2.blast.tsv", "query_length": 30},
    {"path": "passenger-blast/outname_AV_IgYKSTOP_11A_S1_R1_001_atleast-2.blast.tsv", "query_length": 32},
]


def build_processed_df(chigy, blast_paths, query_length):
    blast_df = blast_df_of_blast_files(blast_paths, query_length)
    return chigy.processed_stop_df_of_blast_df(blast_df)


# %% [markdown]
# The pairwise alignment step is expensive (tens of thousands of reads,
# serial Python). We cache the raw `processed_stop_df` outputs to
# parquet in `_ignore/` (already gitignored) so repeated runs of this
# notebook are cheap. Delete `_ignore/mnm_cache_{hc,lc}.parquet` to
# force a rebuild.

# %%
CACHE_DIR = Path("_ignore")
CACHE_DIR.mkdir(exist_ok=True)
HC_CACHE = CACHE_DIR / "mnm_cache_hc.pkl"
LC_CACHE = CACHE_DIR / "mnm_cache_lc.pkl"


# %%
if HC_CACHE.exists():
    print(f"Loading HC from cache {HC_CACHE} ...")
    hc_raw = pd.read_pickle(HC_CACHE)
else:
    print("Building HC processed_stop_df (slow: alignment step) ...")
    hc_raw = build_processed_df(CHIGY_HC, HC_BLAST_PATHS, query_length=20)
    hc_raw.to_pickle(HC_CACHE)
hc_bel = apply_believable_filter(hc_raw)
print(f"HC: {len(hc_raw)} raw, {len(hc_bel)} believable")

# %%
if LC_CACHE.exists():
    print(f"Loading LC from cache {LC_CACHE} ...")
    lc_raw = pd.read_pickle(LC_CACHE)
else:
    print("Building LC processed_stop_df (slow: alignment step) ...")
    lc_frames = []
    for spec in LC_SPECS:
        df = build_processed_df(CHIGY_LC, [spec["path"]], spec["query_length"])
        lc_frames.append(df)
    lc_raw = pd.concat(lc_frames, ignore_index=True)
    lc_raw.to_pickle(LC_CACHE)
lc_bel = apply_believable_filter(lc_raw)
print(f"LC: {len(lc_raw)} raw, {len(lc_bel)} believable")


# %% [markdown]
# ### Assign isotype labels
#
# Runs 5/6/7 have IgM- and IgG-labeled datasets (`_IgM_` / `_IgG_` in
# the dataset name). Runs 8A/9A/11A are IgM-only (they come from the
# `outname_AV_IgYKSTOP_*` series; the "IgY-KSTOP" in the name is the
# passenger construct, not the sequenced isotype), so they are
# pooled into a single "IgM (8A/9A/11A)" panel. Files named
# `outname_[5-7]*_S*_R1_*` without an `_IgM_`/`_IgG_` token are
# dropped from the isotype stratification because their isotype is
# not explicit from the filename.


# %%
def label_isotype(dataset: str) -> str:
    if "_IgM_" in dataset:
        return "IgM"
    if "_IgG_" in dataset:
        return "IgG"
    if "AV_IgYKSTOP" in dataset:
        return "IgM (8A/9A/11A)"
    return "other"


hc_bel["isotype"] = hc_bel["dataset"].map(label_isotype)
lc_bel["isotype"] = lc_bel["dataset"].map(label_isotype)
print("HC isotype counts:")
print(hc_bel["isotype"].value_counts())
print()
print("LC isotype counts:")
print(lc_bel["isotype"].value_counts())


# %% [markdown]
# ## Spatial correlation core
#
# Given a set of reads (each with a sorted list of mutation positions
# on a fixed template of length $L$) and a per-site mutation
# probability vector $p \in [0,1]^L$:
#
# - $n(r)$ — observed pair count at distance $r$: for each read, for
#   each pair $(i < j)$ of mutation positions, increment the bin
#   $r = j - i$.
# - $n_m(r)$ — analytic independent-sites expectation:
#   $n_m(r) = N_\text{reads} \sum_{i<j,\,j-i=r} p_i p_j$.
# - $f(r) = n(r) / n_m(r)$.
#
# Bootstrap CI on $f(r)$: resample reads with replacement `B` times,
# recompute $n(r)$, divide by the (fixed) $n_m(r)$, report
# 2.5 / 97.5 percentiles.
#
# Permutation null band on $f(r)$: for each read with $k$ mutations,
# draw $k$ positions *without replacement* from
# $\{0,\dots,L-1\}$ with probabilities $p_i / \sum_j p_j$, repeat `B`
# times, compute $f(r)$. Under $H_0$ (per-site independence within a
# read), this band should cover 1.

# %%
def observed_pair_histogram(mutation_positions_series, L: int) -> np.ndarray:
    """Return n(r) for r in [0, L-1]; only r >= 1 is meaningful.

    Vectorized: for each read with k mutations, compute all k*(k-1)/2
    pair distances with ``np.subtract.outer``-style broadcasting and
    accumulate into ``hist`` via ``np.bincount``.
    """
    all_diffs = []
    for positions in mutation_positions_series:
        k = len(positions)
        if k < 2:
            continue
        arr = np.sort(np.asarray(positions, dtype=np.int64))
        # j > i so arr[j] - arr[i] > 0
        diffs = arr[None, :] - arr[:, None]
        iu = np.triu_indices(k, k=1)
        all_diffs.append(diffs[iu])
    if not all_diffs:
        return np.zeros(L, dtype=np.int64)
    concat = np.concatenate(all_diffs).astype(np.int64)
    concat = concat[(concat > 0) & (concat < L)]
    hist = np.bincount(concat, minlength=L).astype(np.int64)
    if len(hist) > L:
        hist = hist[:L]
    return hist


def analytic_null_histogram(p: np.ndarray, n_reads: int) -> np.ndarray:
    """n_m(r) = n_reads * sum_{i<j, j-i=r} p_i p_j via FFT-style shifts.

    Implementation: for each shift r in [1, L-1], compute
    sum_i p_i * p_{i+r} and multiply by n_reads.
    """
    L = len(p)
    out = np.zeros(L, dtype=np.float64)
    for r in range(1, L):
        out[r] = float(np.dot(p[:-r], p[r:]))
    return n_reads * out


def empirical_p(bel_df: pd.DataFrame, L: int) -> np.ndarray:
    """Per-site mutation probability (mean over reads).

    Uses ``mutation_frequency_by_position_of`` to get counts/reads,
    then extends with zeros to length L if needed.
    """
    freq = mutation_frequency_by_position_of(bel_df)
    p = np.zeros(L, dtype=np.float64)
    p[: len(freq)] = freq
    return p


def precompute_read_pair_diffs(
    positions_list: list[np.ndarray], L: int
) -> list[np.ndarray]:
    """For each read, return a 1-D int64 array of its pair distances
    (already clipped to ``1 <= d < L``). Reads with <2 mutations give
    an empty array. This lets the bootstrap skip alignment-level work.
    """
    per_read = []
    for positions in positions_list:
        k = len(positions)
        if k < 2:
            per_read.append(np.empty(0, dtype=np.int64))
            continue
        arr = np.sort(np.asarray(positions, dtype=np.int64))
        diffs = arr[None, :] - arr[:, None]
        iu = np.triu_indices(k, k=1)
        d = diffs[iu]
        d = d[(d > 0) & (d < L)]
        per_read.append(d.astype(np.int64))
    return per_read


def hist_from_read_diffs(
    per_read_diffs: list[np.ndarray], indices: np.ndarray, L: int
) -> np.ndarray:
    concat_parts = [per_read_diffs[i] for i in indices]
    if not concat_parts:
        return np.zeros(L, dtype=np.int64)
    concat = np.concatenate(concat_parts) if concat_parts else np.empty(0, dtype=np.int64)
    if concat.size == 0:
        return np.zeros(L, dtype=np.int64)
    return np.bincount(concat, minlength=L)[:L].astype(np.int64)


def bootstrap_f_ci(
    per_read_diffs: list[np.ndarray],
    n_m: np.ndarray,
    L: int,
    n_boot: int = 400,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (f_mean, f_lo, f_hi) from read-level bootstrap using
    precomputed per-read pair-distance arrays."""
    rng = rng or np.random.default_rng()
    n = len(per_read_diffs)
    f_samples = np.zeros((n_boot, L), dtype=np.float64)
    idx_array = np.arange(n)
    for b in range(n_boot):
        idx = rng.choice(idx_array, size=n, replace=True)
        obs = hist_from_read_diffs(per_read_diffs, idx, L)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_samples[b] = np.where(n_m > 0, obs / n_m, np.nan)
    f_mean = np.nanmean(f_samples, axis=0)
    f_lo = np.nanpercentile(f_samples, 2.5, axis=0)
    f_hi = np.nanpercentile(f_samples, 97.5, axis=0)
    return f_mean, f_lo, f_hi


def _gumbel_topk(
    m: int, k: int, log_probs: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Draw ``m`` samples of ``k`` distinct positions from a
    categorical distribution with log-probabilities ``log_probs``
    (length L). Returns an (m, k) array of indices. Uses the Gumbel
    top-k trick: adding Gumbel(0,1) noise to log-probabilities and
    taking the top-k argmax is equivalent to sampling without
    replacement from the softmax distribution.
    """
    L = len(log_probs)
    # Gumbel(0,1): -log(-log(U)) with U ~ Uniform(0,1)
    u = rng.random((m, L))
    # clip to avoid log(0)
    np.clip(u, 1e-12, 1 - 1e-12, out=u)
    gumbel = -np.log(-np.log(u))
    scores = gumbel + log_probs[None, :]
    # top-k indices per row via argpartition
    # argpartition places the k largest in the last k slots (unordered)
    idx = np.argpartition(scores, -k, axis=1)[:, -k:]
    idx.sort(axis=1)
    return idx


def permutation_null_band(
    mutation_counts: np.ndarray,
    p: np.ndarray,
    n_m: np.ndarray,
    L: int,
    n_perm: int = 200,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Null band on f(r) under per-site independence within reads.

    For each read with k mutations, draw k positions without
    replacement from range(L) with probabilities proportional to p,
    build the pair histogram, divide by n_m. Repeat n_perm times.

    Vectorized: groups reads by ``k``, then for each perm and each
    ``k`` draws all M_k reads at once with the Gumbel top-k trick.
    Empty-category reads (k<2) are skipped.
    """
    rng = rng or np.random.default_rng()
    probs = p / p.sum()
    probs = np.clip(probs, 1e-20, None)
    log_probs = np.log(probs)
    counts = mutation_counts.astype(np.int64)
    counts = counts[counts >= 2]
    unique_k, counts_per_k = np.unique(counts, return_counts=True)
    f_samples = np.zeros((n_perm, L), dtype=np.float64)
    for b in range(n_perm):
        all_diffs = []
        for k, m in zip(unique_k, counts_per_k):
            draws = _gumbel_topk(int(m), int(k), log_probs, rng)  # (m, k)
            # pairwise diffs for each row: (m, k, k) - (m, k, k)
            # we only need upper triangle
            diffs = draws[:, None, :] - draws[:, :, None]  # (m, k, k)
            iu = np.triu_indices(int(k), k=1)
            # extract upper-triangle pairs: (m, num_pairs)
            pair_d = diffs[:, iu[0], iu[1]].ravel()
            pair_d = pair_d[(pair_d > 0) & (pair_d < L)]
            all_diffs.append(pair_d.astype(np.int64))
        if all_diffs:
            concat = np.concatenate(all_diffs)
            hist = np.bincount(concat, minlength=L)[:L].astype(np.int64)
        else:
            hist = np.zeros(L, dtype=np.int64)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_samples[b] = np.where(n_m > 0, hist / n_m, np.nan)
    med = np.nanmedian(f_samples, axis=0)
    lo = np.nanpercentile(f_samples, 2.5, axis=0)
    hi = np.nanpercentile(f_samples, 97.5, axis=0)
    return med, lo, hi


# %% [markdown]
# ### Correctness validations (implementation tests)
#
# These must pass before we trust the biological plots. Each test is
# described in the issue's "Correctness validations" section.

# %% [markdown]
# **Test 3 — two-mutation degenerate case.** For reads with exactly 2
# mutations, the total pair count equals the read count.

# %%
def test_two_mutation_case():
    mps = [np.array([3, 17]), np.array([1, 2]), np.array([5, 29])]
    hist = observed_pair_histogram(mps, L=100)
    assert hist.sum() == 3, (hist.sum(), 3)
    # expected pair distances: 14, 1, 24
    assert hist[14] == 1
    assert hist[1] == 1
    assert hist[24] == 1
    print("test_two_mutation_case: OK")


test_two_mutation_case()


# %% [markdown]
# **Test 1 — known-answer synthetic null.** The analytic null formula
# $n_m(r) = N \sum p_i p_j$ assumes *per-site independent Bernoulli*
# draws (each site `i` mutates independently with probability `p_i`).
# To validate the analytic formula end-to-end, we generate synthetic
# reads under exactly that model and check that $f(r) \approx 1$.
#
# Note that this is **not** the same as drawing a fixed `k` positions
# without replacement — the two differ for concentrated `p` (e.g. HC
# hotspots). The per-site Bernoulli generator gives variable `k` per
# read but exactly matches the analytic null's implicit generative
# model. The "fixed-k" perturbation on real data is handled by the
# permutation null band in the main analysis below, which uses the
# per-read `k` distribution as a constraint.

# %%
def generate_bernoulli_reads(
    n_reads: int, p: np.ndarray, rng: np.random.Generator
) -> list[np.ndarray]:
    """Draw ``n_reads`` synthetic reads under per-site Bernoulli(p_i)."""
    L = len(p)
    # draw a big uniform matrix and compare to p
    uniforms = rng.random((n_reads, L))
    mask = uniforms < p[None, :]
    out = []
    for row in mask:
        out.append(np.where(row)[0].astype(np.int64))
    return out


def generate_synthetic_reads_fixed_k(
    mutation_counts: np.ndarray,
    p: np.ndarray,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Draw ``len(mutation_counts)`` reads: each has `k` positions
    sampled without replacement with probabilities proportional to
    `p`. Used by the permutation null band, NOT by the analytic-null
    correctness test."""
    L = len(p)
    probs = p / p.sum()
    positions = np.arange(L)
    out = []
    for k in mutation_counts:
        k = int(k)
        if k < 2:
            out.append(np.array([], dtype=np.int64))
            continue
        draw = rng.choice(positions, size=k, replace=False, p=probs)
        draw.sort()
        out.append(draw)
    return out


def run_known_answer_null(
    bel_df: pd.DataFrame, L: int, label: str, n_synth: int = 50000
) -> None:
    p = empirical_p(bel_df, L)
    synth_positions = generate_bernoulli_reads(n_synth, p, RNG)
    obs = observed_pair_histogram(synth_positions, L)
    n_m = analytic_null_histogram(p, n_reads=n_synth)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(n_m > 0, obs / n_m, np.nan)
    r_vals = np.arange(1, min(31, L))
    f_slice = f[r_vals]
    mean_f = float(np.nanmean(f_slice))
    max_dev = float(np.nanmax(np.abs(f_slice - 1)))
    print(f"{label} synthetic null: mean f(r) over r∈[1,30] = {mean_f:.3f}")
    print(f"  min {np.nanmin(f_slice):.3f}, max {np.nanmax(f_slice):.3f}")
    print(f"  max |f(r) - 1| = {max_dev:.3f}")
    # tolerance: with 50k synthetic reads, Poisson noise on n(r) is
    # ~O(1/sqrt(n_m(r))) which is small; a passing test should be
    # within ~10% of 1 at short r.
    assert mean_f < 1.1 and mean_f > 0.9, (
        f"{label} synthetic null mean f(r)={mean_f:.3f} far from 1; "
        "the analytic null formula or the histogram is buggy"
    )


run_known_answer_null(hc_bel, L=len(CHIGY_HC.chigy_stop_trimmed), label="HC")
run_known_answer_null(lc_bel, L=len(CHIGY_LC.chigy_stop_trimmed), label="LC")


# %% [markdown]
# **Test 2 — within-read position shuffle invariance.** For each real
# read, replace its positions by $k$ positions drawn without
# replacement from the empirical marginal. This gives a null in which
# co-localization is destroyed but the per-read `k` distribution is
# preserved. Computed against the analytic null, $f(r)$ should be
# approximately constant in $r$ (no short-range excess), though the
# absolute level will differ from 1 due to the fixed-`k` constraint
# (the analytic null assumes Bernoulli per site). The point of this
# test is that the short-range *shape* collapses — if real data shows
# `f(r)` decreasing in `r` and the shuffle shows no trend, the signal
# is structural, not an artifact of the histogram code.

# %%
def run_shuffle_invariance(bel_df: pd.DataFrame, L: int, label: str) -> None:
    counts = bel_df["mutation_count"].to_numpy()
    p = empirical_p(bel_df, L)
    shuffled = generate_synthetic_reads_fixed_k(counts, p, RNG)
    obs = observed_pair_histogram(shuffled, L)
    n_m = analytic_null_histogram(p, n_reads=len(counts))
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(n_m > 0, obs / n_m, np.nan)
    r_vals = np.arange(1, min(31, L))
    f_slice = f[r_vals]
    # measure "flatness" as max - min over the range
    flatness = float(np.nanmax(f_slice) - np.nanmin(f_slice))
    mean_f = float(np.nanmean(f_slice))
    print(
        f"{label} shuffle: mean f(r) = {mean_f:.3f}, "
        f"range {float(np.nanmin(f_slice)):.3f} – {float(np.nanmax(f_slice)):.3f}, "
        f"flatness = {flatness:.3f}"
    )


run_shuffle_invariance(hc_bel, L=len(CHIGY_HC.chigy_stop_trimmed), label="HC")
run_shuffle_invariance(lc_bel, L=len(CHIGY_LC.chigy_stop_trimmed), label="LC")


# %% [markdown]
# ## Main plots: $f(r)$ on real data
#
# For each (chain, stratum) pair, compute $f(r)$ on the real believable
# reads, overlay the bootstrap 95% CI and the permutation null 95%
# band. Stratify the HC panel by IgM vs IgG (runs 5/6/7 only, plus a
# pooled "all-HC" panel). LC is single-isotype (IgM, pooled over
# 8A/9A/11A).

# %%
R_MAX = 30


def compute_fr_stats(
    bel_df: pd.DataFrame, L: int, n_boot: int = 400, n_perm: int = 200
):
    if len(bel_df) == 0:
        raise ValueError("empty dataframe")
    p = empirical_p(bel_df, L)
    positions_list = [np.asarray(x, dtype=np.int64) for x in bel_df["mutation_positions"]]
    counts = bel_df["mutation_count"].to_numpy()
    n_m = analytic_null_histogram(p, n_reads=len(bel_df))
    obs = observed_pair_histogram(positions_list, L)
    with np.errstate(divide="ignore", invalid="ignore"):
        f_point = np.where(n_m > 0, obs / n_m, np.nan)
    per_read_diffs = precompute_read_pair_diffs(positions_list, L)
    f_mean, f_lo, f_hi = bootstrap_f_ci(per_read_diffs, n_m, L, n_boot=n_boot, rng=RNG)
    null_med, null_lo, null_hi = permutation_null_band(
        counts, p, n_m, L, n_perm=n_perm, rng=RNG
    )
    return dict(
        f_point=f_point,
        f_lo=f_lo,
        f_hi=f_hi,
        null_lo=null_lo,
        null_hi=null_hi,
        null_med=null_med,
        n_reads=len(bel_df),
    )


def plot_fr(ax, stats, label: str, r_max: int = R_MAX):
    r = np.arange(1, r_max + 1)
    ax.axhline(1, color="k", lw=0.5, ls="--")
    ax.fill_between(
        r, stats["null_lo"][r], stats["null_hi"][r], alpha=0.2, color="gray",
        label="perm null 95%",
    )
    ax.fill_between(
        r, stats["f_lo"][r], stats["f_hi"][r], alpha=0.3, color="C0",
        label="bootstrap 95%",
    )
    ax.plot(r, stats["f_point"][r], color="C0", lw=1.5)
    ax.set_xlabel("pair distance r (nt)")
    ax.set_ylabel("f(r) = n(r) / n_m(r)")
    ax.set_title(f"{label}  (N={stats['n_reads']})")


# %% [markdown]
# ### HC: all runs pooled, stratified by isotype and mutation load

# %%
def hc_panels():
    L = len(CHIGY_HC.chigy_stop_trimmed)
    panels = {}
    panels["HC all, ≤9 muts"] = hc_bel
    panels["HC all, ≤5 muts"] = hc_bel[hc_bel["mutation_count"] <= LOW_LOAD_CAP]
    for iso in ["IgM", "IgG"]:
        sub = hc_bel[hc_bel["isotype"] == iso]
        if len(sub) > 0:
            panels[f"HC {iso}, ≤9 muts"] = sub
            panels[f"HC {iso}, ≤5 muts"] = sub[sub["mutation_count"] <= LOW_LOAD_CAP]
    stats = {k: compute_fr_stats(v, L) for k, v in panels.items()}
    return stats


print("Computing HC panels ...")
hc_stats = hc_panels()

# %%
fig, axes = plt.subplots(
    nrows=(len(hc_stats) + 1) // 2, ncols=2, figsize=(12, 3 * ((len(hc_stats) + 1) // 2)),
    sharex=True,
)
axes = np.asarray(axes).flatten()
for ax, (label, stats) in zip(axes, hc_stats.items()):
    plot_fr(ax, stats, label)
    ax.set_ylim(0, max(4, np.nanmax(stats["f_hi"][1 : R_MAX + 1]) * 1.1))
for ax in axes[len(hc_stats) :]:
    ax.axis("off")
axes[0].legend(loc="upper right", fontsize=8)
fig.suptitle("Heavy chain passenger  —  spatial correlation f(r)")
fig.tight_layout()


# %% [markdown]
# ### LC: pooled IgM (8A / 9A / 11A)

# %%
def lc_panels():
    L = len(CHIGY_LC.chigy_stop_trimmed)
    panels = {
        "LC pooled, ≤9 muts": lc_bel,
        "LC pooled, ≤5 muts": lc_bel[lc_bel["mutation_count"] <= LOW_LOAD_CAP],
    }
    stats = {k: compute_fr_stats(v, L) for k, v in panels.items()}
    return stats


print("Computing LC panels ...")
lc_stats = lc_panels()

# %%
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), sharex=True)
for ax, (label, stats) in zip(axes, lc_stats.items()):
    plot_fr(ax, stats, label)
    ax.set_ylim(0, max(4, np.nanmax(stats["f_hi"][1 : R_MAX + 1]) * 1.1))
axes[0].legend(loc="upper right", fontsize=8)
fig.suptitle("Light chain passenger  —  spatial correlation f(r)")
fig.tight_layout()


# %% [markdown]
# ## Success-criteria readout
#
# Declare **positive** if, in at least one chain, the bootstrap lower
# bound of $f(r)$ exceeds 1.5 at some $r \in [2, 5]$.

# %%
def summarize_panel(label: str, stats: dict) -> dict:
    f_point = stats["f_point"]
    f_lo = stats["f_lo"]
    f_hi = stats["f_hi"]
    r_short = np.arange(2, 6)
    lo_short = f_lo[r_short]
    point_short = f_point[r_short]
    return {
        "panel": label,
        "n_reads": stats["n_reads"],
        "f_mean_r2_5": float(np.nanmean(point_short)),
        "f_lo_min_r2_5": float(np.nanmin(lo_short)),
        "f_lo_max_r2_5": float(np.nanmax(lo_short)),
        "positive_at_r2_5": bool(np.nanmax(lo_short) > 1.5),
    }


summary_rows = []
for label, stats in {**hc_stats, **lc_stats}.items():
    summary_rows.append(summarize_panel(label, stats))
summary_df = pd.DataFrame(summary_rows)
summary_df

# %%
any_positive = summary_df["positive_at_r2_5"].any()
print("=" * 60)
if any_positive:
    print("RESULT: POSITIVE — at least one panel shows f(r) bootstrap")
    print("lower bound > 1.5 at some r ∈ [2, 5].")
else:
    print("RESULT: NEGATIVE or AMBIGUOUS — no panel crosses the")
    print("positive threshold. Inspect plots and reconsider caps /")
    print("stratification.")
print("=" * 60)
