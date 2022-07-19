import pickle, argparse, os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "--treedir", type=str, help="directory with pickled gctrees to read"
)
parser.add_argument("--out", type=str, help="csv outfile")
parser.add_argument(
    "--obs",
    action="store_true",
    help="display statistics only for observed (leaf) nodes",
)
args = parser.parse_args()


def replay_tree_summary_stats(tree):
    sequences = []
    distances = []
    abundances = []
    isotypes = []
    LBIs = []
    affinities = []
    if args.obs:
        nodes = tree.get_leaves()
    else:
        nodes = list(tree.traverse())
    for node in nodes:
        sequences.append(node.sequence)
        distances.append(node.get_distance(tree))
        abundances.append(node.abundance)
        isotypes.extend(isotype.iso_string() for isotype in node.isotype.keys())
        LBIs.append(node.LBI)
        affinities.append(node.delta_bind_CGG_tdms_linear_model_pred)
    clade_sizes = [
        sum(node.abundance for node in child.traverse()) for child in tree.children
    ]
    norm_dominance_score = max(clade_sizes) / sum(clade_sizes)
    return (
        len(sequences),
        np.median(abundances),
        np.mean(abundances),
        np.mean(distances),
        np.median(LBIs),
        np.mean(LBIs),
        np.min(LBIs),
        np.max(LBIs),
        norm_dominance_score,
        np.median(affinities),
        np.mean(affinities),
        np.min(affinities),
        np.max(affinities),
        np.std(affinities),
        100 * isotypes.count("IgG1") / len(isotypes),
        100 * isotypes.count("IgG2") / len(isotypes),
        100 * isotypes.count("IgG3") / len(isotypes),
        100 * isotypes.count("IgM") / len(isotypes),
        100 * isotypes.count("IgD") / len(isotypes),
        100 * isotypes.count("IgA") / len(isotypes),
        100 * isotypes.count("IgE") / len(isotypes),
    )


stats_dict = {}

for file in os.listdir(args.treedir):
    if file.endswith("p"):
        with open(f"{args.treedir}/{file}", "rb") as fh:
            tree = pickle.load(fh)
            stats_dict[file] = replay_tree_summary_stats(tree.tree)

# write csv
with open(args.out, "w") as fh:
    print(
        f"tree filename,number of unique sequences,median abundance,mean abundance,mean number of substitutions,median LBI,mean LBI,min LBI,max LBI,normalized dominance score, median affinity, mean affinity,min affinity,max affinity, affinity std, %IgG1,%IgG2,%IgG3,%IgM,%IgD,%IgA,%IgE",
        file=fh,
    )
    for name, stats in stats_dict.items():
        print(f'{name},{",".join(str(stat) for stat in stats)}', file=fh)
