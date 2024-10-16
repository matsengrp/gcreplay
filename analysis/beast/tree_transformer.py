from lark import Transformer, Lark
from lark.visitors import Discard
import ete3


class TreelistTransformer(Transformer):
    list = list

    def integer(self, s):
        (s,) = s
        return int(s)

    def float(self, s):
        (s,) = s
        return float(s)

    def string(self, s):
        (s,) = s
        return str(s)

    def support(self, s):
        return [("dist", s[0])]

    def quoted_escaped_string(self, s):
        (value,) = s
        return str(value[1:-1])

    def annotationname(self, s):
        return str(s[0])

    def annotation(self, s):
        return tuple(s)

    def annotationlist(self, s):
        return list(s)

    def extendedannotationlist(self, s):
        return [it for ls in s for it in ls]

    def internal_node(self, s):
        return self.leaf_node(s)

    def leaf_node(self, s):
        n = ete3.TreeNode()
        for parse_child in s:
            if parse_child.data == "children":
                for child in parse_child.children:
                    n.add_child(child=child)
            elif parse_child.data == "annotations":
                n.add_features(
                    **dict(t for child in parse_child.children for t in child)
                )
            elif parse_child.data == "nodename":
                n.name = str(parse_child.children[0])
        return n


class NexusTransformer(TreelistTransformer):
    def pair(self, s):
        return tuple(s)

    def taxon_name(self, s):
        return str(s[0])

    def leaf_mapping(self, s):
        leaf_map = dict(s)
        self.leaf_map = leaf_map
        return leaf_map

    def translate_block_wrapper(self, s):
        return Discard

    def leaf_node(self, s):
        n = ete3.TreeNode()
        for parse_child in s:
            if parse_child.data == "children":
                for child in parse_child.children:
                    n.add_child(child=child)
            elif parse_child.data == "annotations":
                n.add_features(
                    **dict(t for child in parse_child.children for t in child)
                )
            elif parse_child.data == "nodename":
                n.name = self.leaf_map[int(parse_child.children[0])]
        return n

    def beast_tree(self, s):
        tree = s[-1]
        tree.add_feature("beast_name", s[0])
        tree.add_feature("beast_annotations", dict(s[1]))
        tree.add_feature("rice_annotation", s[2])
        return tree

    def beast_file(self, s):
        return s[1]


with open("newick.lark", "r") as fh:
    newick_lark_string = fh.read()

_nexus_parser = Lark(
    newick_lark_string,
    parser="lalr",
    start="beast_file",
    transformer=NexusTransformer(),
)

_nexus_translate_block_parser = Lark(
    newick_lark_string,
    parser="lalr",
    start="translate_block",
    transformer=NexusTransformer(),
)

_newick_parser = Lark(
    newick_lark_string,
    parser="lalr",
    start="newicklist",
    transformer=TreelistTransformer(),
)


class NexusIterator:
    def __init__(self, nexus_filename):
        self.file = nexus_filename
        self.nexus_newick_parser = self.get_beast_tree_parser()

    def get_beast_tree_parser(self):
        with open(self.file, "r") as fh:
            for line in fh:
                if "Begin trees;" in line:
                    break
            data = ""
            for line in fh:
                if ";" in line:
                    data += line
                    break
                else:
                    data += line

        transformer = NexusTransformer()
        transformer.leaf_map = _nexus_translate_block_parser.parse(data)
        return Lark(
            newick_lark_string,
            parser="lalr",
            start="beast_tree",
            transformer=transformer,
        )

    def iter_tree_sections(self):
        with open(self.file, "r") as fh:
            for line in fh:
                if "tree STATE" in line[:12]:
                    current_tree = line
                    break
            for line in fh:
                if "tree STATE" in line[:12]:
                    yield current_tree
                    current_tree = line
                elif "End;" in line[:12]:
                    break
                else:
                    current_tree += line
        yield current_tree

    def iter_trees(self):
        for treestring in self.iter_tree_sections():
            yield self.nexus_newick_parser.parse(treestring)


def parse_newick(data):
    """Parse a string containing concatenated newick strings.

    Returns a list of ete3.Tree objects containing the parsed trees."""
    return _newick_parser.parse(data)


def parse_newick_file(filepath):
    """Parse a file containing newick strings.

    Returns a list of ete3.Tree objects containing the parsed trees."""
    with open(filepath, "r") as fh:
        return _newick_parser.parse(fh.read())


def parse_nexus_file(filepath):
    """Parse a nexus file output by Beast.

    Returns a list of ete3.Tree objects containing the parsed trees. Leaves
    have taxon names in the ``name`` attribute, and tree root nodes have
    ``beast_name``, ``beast_annotations``, and ``rice_annotation``
    attributes containing tree annotations. All other extended newick
    annotations are saved in the corresponding attribute on each node.
    """
    with open(filepath, "r") as fh:
        for line in fh:
            if "Begin trees;" in line:
                data = fh.read()
    return _nexus_parser.parse(data)


def iter_nexus_trees(filepath):
    """Parse a nexus file output by Beast.
    Loads and yields one tree at a time, avoiding loading all trees in memory
    at once.
    """
    parser_wrapper = NexusIterator(filepath)
    yield from parser_wrapper.iter_trees()
