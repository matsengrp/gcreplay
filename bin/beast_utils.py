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
        partial(beast_dendropy_to_ete, phenotype_fn=phenotype_fn), 
        beast_trees
    )
    bad_root = lambda tree: not tree.get_leaves_by_name(naive_label)[0].up.is_root()
    return [tree for tree in ete_iterator if not bad_root(tree)]


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
    tree_start_idx = int(len(beast_trees)*burn_frac)
    any((set_dendropy_tree_node_ages(tree) for tree in beast_trees[tree_start_idx:]))
    return tree_start_idx, beast_trees[tree_start_idx:]


#def set_dendropy_tree_node_ages(tree):
#    """
#    Set the age and root_distance attributes (in place) of nodes in tree. The
#    age of a node is taken to be longest path from the root to any leaf minus
#    the length of the path from the root to the given node. The age of non-naive
#    leaf nodes may not be exactly 0.0 due to floating point arithmetic.
#    """
#    last_time = max(tree.calc_node_root_distances())
#    def age_function(node): 
#        age = last_time - node.root_distance
#        if node.is_leaf():
#            age = round(age, 6)
#            assert age < node.up.age
#        return age
#    tree.calc_node_ages(set_node_age_fn=age_function)
#    return None


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
        # if len(dendropy_node.comments) == 2:
        #     mutation_history = format_history_all(dendropy_node.comments[1])
        #     mutation_history.sort(key=lambda x: x[1], reverse=True)
        #     parent_sequence = list(ete_node.up.sequence)

        #     previous_event_age = ete_node.up.age
        #     for site, age, ancestral, derived in mutation_history:
        #         site -= 1
        #         time_to_event = previous_event_age - age
        #         parent_sequence[site] = derived

        #         event = dict(
        #                 site=site,
        #                 ancestral=ancestral,
        #                 derived=derived,
        #                 age=age,
        #                 time_to_event=time_to_event,
        #                 phenotypes=phenotype_fn("".join(parent_sequence))
        #         )

        #         ete_node.mutations.append(event)
        #         previous_event_age = age

    # now, we want to fix the tree such that naive is the root
    # outgroup = tree.children[1]
    outgroup = ete_tree&"naive@0"
    ete_tree.remove_child(outgroup)
    ete_tree.dist = outgroup.dist
    outgroup.dist=0
    outgroup.add_child(ete_tree)
    ete_tree = outgroup

    for node in ete_tree.traverse():
        node.time = node.get_distance(ete_tree)

    leaf_ages = [leaf.time for leaf in ete_tree]
    assert np.allclose(leaf_ages[0], leaf_ages)
    time_scale_factor = leaf_ages[0]
    for node in ete_tree.traverse():
        node.time /= time_scale_factor
        node.dist /= time_scale_factor
        del node.age

    assert len(ete_tree.children) == 1
    # TODO check that I should be setting these below
    # ete_tree.children[0].dist = 0
    # ete_tree.children[0].time = 0
    assert ete_tree.children[0].time == 0, f"naive time = {ete_tree.children[0].time}"
    for node in ete_tree.traverse():
        assert np.isclose(node.get_distance(ete_tree), node.time), f"node {node.name} has distance {node.get_distance(ete_tree)} != time {node.time}"
        if node.is_leaf():
            assert node.time == 1, f"leaf time should be 1, {node.time}"    

    return ete_tree
