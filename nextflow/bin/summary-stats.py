import pickle, argparse, os
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('--treedir',type=str)
parser.add_argument('--out',type=str)
args = parser.parse_args()


def replay_tree_summary_stats(tree):
    sequences = []
    distances = []
    abundances = []
    isotypes = []
    LBIs = []
    affinities = []

    for node in tree.traverse():
        sequences.append(node.sequence)
        distances.append(node.get_distance(tree))
        abundances.append(node.abundance)
        isotypes.append((isotype.iso_string() for isotype in node.isotype.keys()))
        LBIs.append(node.LBI)
        affinities.append(node.delta_bind_CGG_tdms_linear_model_pred)
    clade_sizes = [sum(node.abundance for node in child.traverse()) for child in tree.children]
    norm_dominance_score = max(clade_sizes) / sum(clade_sizes)
    return len(sequences), np.median(abundances), np.mean(abundances), np.mean(distances), np.median(LBIs), np.mean(LBIs), np.min(LBIs), np.max(LBIs), norm_dominance_score, np.median(affinities), np.mean(affinities), np.min(affinities), np.max(affinities), np.std(affinities)

stats_dict = {}

for file in os.listdir(args.treedir):
    if file.endswith("p"):
        with open(f"{args.treedir}/{file}", 'rb') as fh:
            tree = pickle.load(fh)
            stats_dict[file] = replay_tree_summary_stats(tree.tree)

# write csv
with open(args.out, 'w') as fh:
    print(f'tree filename,number of unique sequences,median abundance,mean abundance,mean number of substitutions,median LBI,mean LBI,min LBI,max LBI,normalized dominance score, median affinity, mean affinity,min affinity,max affinity, affinity std', file=fh)
    for name, stats in stats_dict.items():
        print(f'{name},{",".join(str(stat) for stat in stats)}', file=fh)
