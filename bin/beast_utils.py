import ete3
from historydag import beast_loader
from functools import partial
import numpy as np


def beast_trees_as_ete(
    beast_xml_path,
    beast_nexus_path,
    naive_sequence,
    phenotype_fn,
    naive_label="naive@0",
    skip_first_tree=True,
    burn_frac=0,
):
    """
    Returns a list of ete3 Tree objects from BEAST files. The trees are fully
    processed so that nodes of trees have attributes age, mutations, and
    sequence required by complex_subs_by_fivemer.log_probability_for_ete_tree.
    Trees that do not have the naive seqeunce as an outgroup are omitted. By
    default, the first tree of the beast file (tree state 0) is omitted.
    """
    tsi, beast_trees = dendropy_trees_from_beast(
        beast_xml_path, beast_nexus_path, naive_sequence, skip_first_tree, burn_frac
    )
    return tsi, dendropy_beast_trees_to_ete(beast_trees, naive_label, phenotype_fn)


def dendropy_beast_trees_to_ete(beast_trees, naive_label, phenotype_fn):
    """
    Returns a list of ete3 Tree objects from a list of dendropy tree objects.
    The dendropy trees are as described in the documentation for
    dendropy_trees_from_beast; the ete3 trees are as described in the
    documentation for beast_trees_as_ete.
    """
    ete_iterator = map(
        partial(beast_dendropy_to_ete, phenotype_fn=phenotype_fn), beast_trees
    )
    return [tree for tree in ete_iterator]


def dendropy_trees_from_beast(
    beast_xml_path, beast_nexus_path, naive_sequence, skip_first_tree, burn_frac
):
    """
    Returns a list of dendropy Tree objects from BEAST files, which are
    loaded and process by the beast loaded class of historydag. The nodes
    of trees have attributes age and comments, the latter of which contains
    sequence and mutation history.
    """
    tree_iterator = beast_loader.load_beast_trees(
        beast_xml_path, beast_nexus_path, reference_sequence=naive_sequence
    )[0]
    if skip_first_tree:
        next(tree_iterator)
    beast_trees = list(tree_iterator)
    tree_start_idx = int(len(beast_trees) * burn_frac)
    any((set_dendropy_tree_node_ages(tree) for tree in beast_trees[tree_start_idx:]))
    return tree_start_idx, beast_trees[tree_start_idx:]


def set_dendropy_tree_node_ages(tree):
    """
    Set the age and root_distance attributes (in place) of nodes in tree. The
    age of a node is taken to be longest path from the root to any leaf minus
    the length of the path from the root to the given node. The age of non-naive
    leaf nodes may not be exactly 0.0 due to floating point arithmetic.
    """
    last_time = max(tree.calc_node_root_distances())
    age_function = lambda node: last_time - node.root_distance
    tree.calc_node_ages(set_node_age_fn=age_function)

    for node in tree.leaf_node_iter():
        node.age = round(node.age)

    return None


def format_history_all(history_all_string):
    """
    Process the '&history_all' string from BEAST, which should be the second
    entry of a dendropy tree node's comments attribute, and return a list for
    the mutations attribute of the corresponding node in an ete tree.
    """
    return [
        [int((x := mutation[1:].split(","))[0]), float(x[1]), *x[2:]]
        for mutation in history_all_string[14:-2].split("},")
    ]


def beast_dendropy_to_ete(tree, phenotype_fn):
    """
    Return an ete3 tree, whose nodes have age, mutations, and sequence
    attributes, based on the dendropy tree.
    """

    # assert False

    assert tree.is_rooted
    ete_tree = ete3.Tree(dist=0.0)
    node_pairing = zip(
        tree.preorder_node_iter(), ete_tree.traverse(strategy="preorder")
    )
    for dendropy_node, ete_node in node_pairing:
        ete_node.sequence = dendropy_node.comments[0].strip('&states="')
        ete_node.age = dendropy_node.age
        ete_node.mutations = []
        ete_node.phenotypes = phenotype_fn(ete_node.sequence)

        for child in dendropy_node.child_nodes():
            node_label = child.taxon.label if child.taxon else None
            ete_child = ete3.Tree(name=node_label, dist=child.edge_length)
            ete_node.add_child(ete_child)

        # assert False
        if len(dendropy_node.comments) == 2:
            # assert False, f"comments[1]: \n {dendropy_node.comments[1]}"
            mutation_history = format_history_all(dendropy_node.comments[1])
            mutation_history.sort(key=lambda x: x[1], reverse=True)
            parent_sequence = list(ete_node.up.sequence)

            # previous_event_age = ete_node.up.age
            for site, age, ancestral, derived in mutation_history:
                site -= 1
                # time_to_event = previous_event_age - age
                time_since_parent = ete_node.up.age - age
                time_until_child = age - ete_node.age

                ancestral_phenotypes = phenotype_fn("".join(parent_sequence))

                parent_sequence[site] = derived

                event = dict(
                    site=site,
                    ancestral=ancestral,
                    derived=derived,
                    # age=age,
                    time_since_parent=time_since_parent,
                    time_until_child=time_until_child,
                    derived_phenotypes=phenotype_fn("".join(parent_sequence)),
                    ancestral_phenotypes=ancestral_phenotypes,
                )

                ete_node.mutations.append(event)
                # previous_event_age = age

    # now, we want to fix the tree such that naive is the root
    outgroup = ete_tree & "naive@0"

    # make sure the psuedo root has no mutations
    assert not ete_tree.mutations

    # re-root the tree such that naive is new root.
    ete_tree.remove_child(outgroup)
    ete_tree.dist = outgroup.dist
    outgroup.dist = 0
    outgroup.add_child(ete_tree)
    ete_tree = outgroup

    # now, we need to correctly add mutations that were from pseudo root to naive
    # as mutations from naive to pseudo root
    psuedo_root = ete_tree.get_children()[0]
    if ete_tree.mutations:
        for mut_event in reversed(ete_tree.mutations):
            reversion_event = dict(
                site=mut_event["site"],
                ancestral=mut_event["derived"],
                derived=mut_event["ancestral"],
                time_since_parent=mut_event["time_until_child"],
                time_until_child=mut_event["time_since_parent"],
                derived_phenotypes=mut_event["ancestral_phenotypes"],
                ancestral_phenotypes=mut_event["derived_phenotypes"],
            )

            psuedo_root.mutations.append(reversion_event)

    # now remove the mutations from the new naive root.
    ete_tree.mutations = []

    # rescale times to be between 0 and 1
    for node in ete_tree.traverse():
        node.time = node.get_distance(ete_tree)

    leaf_ages = [leaf.time for leaf in ete_tree]
    assert np.allclose(leaf_ages[0], leaf_ages)
    time_scale_factor = leaf_ages[0]
    for node in ete_tree.traverse(strategy="preorder"):
        node.time /= time_scale_factor
        node.dist /= time_scale_factor

        del node.age

        if node.mutations:
            for mut_event in node.mutations:
                mut_event["time_since_parent"] /= time_scale_factor
                mut_event["time_until_child"] /= time_scale_factor
                mut_event["time"] = node.up.time + mut_event["time_since_parent"]
                assert np.isclose(
                    mut_event["time"], node.time - mut_event["time_until_child"]
                )

                # we no longer need these after correcting root mutations
                del mut_event["time_since_parent"]
                del mut_event["time_until_child"]

                assert mut_event["time"] > node.up.time
                assert mut_event["time"] < node.time

    # some final checks on the tree
    assert len(ete_tree.children) == 1
    assert ete_tree.time == 0, f"naive time = {ete_tree.time}"
    for node in ete_tree.traverse():
        assert np.isclose(
            node.get_distance(ete_tree), node.time
        ), f"node {node.name} has distance {node.get_distance(ete_tree)} != time {node.time}"
        if node.is_leaf():
            assert np.isclose(node.time, 1), f"leaf time should be 1, {node.time}"

        for i in range(len(node.mutations) - 1):
            assert (
                node.mutations[i]["time"] < node.mutations[i + 1]["time"]
            ), f"mutation {i} time {node.mutations[i]['time']} should be less than mutation {i+1} time {node.mutations[i+1]['time']}"

    return ete_tree
