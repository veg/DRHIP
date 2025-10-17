"""
Tree manipulation utilities for HyPhy results.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import re
from collections import Counter
from typing import Dict, Optional

# Codon translation table
TT_table = {
    "TTT": "F",
    "TCT": "S",
    "TAT": "Y",
    "TGT": "C",
    "TTC": "F",
    "TCC": "S",
    "TAC": "Y",
    "TGC": "C",
    "TTA": "L",
    "TCA": "S",
    "TAA": "*",
    "TGA": "*",
    "TTG": "L",
    "TCG": "S",
    "TAG": "*",
    "TGG": "W",
    "CTT": "L",
    "CCT": "P",
    "CAT": "H",
    "CGT": "R",
    "CTC": "L",
    "CCC": "P",
    "CAC": "H",
    "CGC": "R",
    "CTA": "L",
    "CCA": "P",
    "CAA": "Q",
    "CGA": "R",
    "CTG": "L",
    "CCG": "P",
    "CAG": "Q",
    "CGG": "R",
    "ATT": "I",
    "ACT": "T",
    "AAT": "N",
    "AGT": "S",
    "ATC": "I",
    "ACC": "T",
    "AAC": "N",
    "AGC": "S",
    "ATA": "I",
    "ACA": "T",
    "AAA": "K",
    "AGA": "R",
    "ATG": "M",
    "ACG": "T",
    "AAG": "K",
    "AGG": "R",
    "GTT": "V",
    "GCT": "A",
    "GAT": "D",
    "GGT": "G",
    "GTC": "V",
    "GCC": "A",
    "GAC": "D",
    "GGC": "G",
    "GTA": "V",
    "GCA": "A",
    "GAA": "E",
    "GGA": "G",
    "GTG": "V",
    "GCG": "A",
    "GAG": "E",
    "GGG": "G",
    "---": "-",
    "NNN": "?",
}
nucs = set(["A", "C", "G", "T"])


def TT(codon: str) -> str:
    """Translate a codon to an amino acid.

    Args:
        codon: A three-letter nucleotide codon

    Returns:
        The corresponding amino acid or '?' for invalid codons
    """
    if codon in TT_table:
        return TT_table[codon]
    return "?"


def newick_parser(
    nwk_str: str,
    bootstrap_values: bool,
    track_tags: Optional[Dict] = None,
    optional_starting_tags: Optional[Dict] = None,
) -> Dict:
    """Parse a Newick format tree string into a hierarchical JSON structure.

    Args:
        nwk_str: Newick format tree string
        bootstrap_values: Whether to interpret node names as bootstrap values for internal nodes
        track_tags: Optional dictionary to track node tags
        optional_starting_tags: Optional dictionary of starting tags

    Returns:
        Dictionary containing the parsed tree or an error message
    """
    if optional_starting_tags is None:
        optional_starting_tags = {}

    clade_stack = []
    automaton_state = 0
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
    quote_delimiter = None
    name_quotes = {"'": 1, '"': 1}

    def add_new_tree_level():
        new_level = {"name": None}
        the_parent = clade_stack[len(clade_stack) - 1]
        if "children" not in the_parent:
            the_parent["children"] = []

        clade_stack.append(new_level)
        the_parent["children"].append(clade_stack[len(clade_stack) - 1])
        clade_stack[len(clade_stack) - 1]["original_child_order"] = len(
            the_parent["children"]
        )

    def finish_node_definition():
        nonlocal current_node_name
        nonlocal current_node_annotation
        nonlocal current_node_attribute

        this_node = clade_stack.pop()
        if bootstrap_values:
            # For the root node, set bootstrap_values at root level
            if len(clade_stack) == 0:  # This is the root node
                tree_json["bootstrap_values"] = current_node_name
            elif "children" in this_node:  # For internal nodes
                this_node["bootstrap_values"] = current_node_name
            else:  # For leaf nodes
                this_node["name"] = current_node_name
        else:
            this_node["name"] = current_node_name

        this_node["attribute"] = current_node_attribute
        this_node["annotation"] = current_node_annotation

        try:
            if "children" not in this_node:
                node_tag = "background"
            for k, v in optional_starting_tags.items():
                if this_node["name"].find(k) >= 0:
                    node_tag = v
                    break
            else:
                node_tag = "test"

            this_node["tag"] = node_tag
        except Exception as exc:
            print(f"Exception {exc}")

        if track_tags is not None:
            track_tags[this_node["name"]] = this_node["tag"]

        current_node_name = ""
        current_node_attribute = ""
        current_node_annotation = ""

    def generate_error(location):
        return {
            "json": None,
            "error": f"Unexpected '{nwk_str[location]}' in '{nwk_str[location - 20 : location + 1]}[ERROR HERE]{nwk_str[location + 1 : location + 20]}'",
        }

    tree_json = {"name": "root"}

    clade_stack.append(tree_json)

    space = re.compile(r"\s")

    for char_index in range(len(nwk_str)):
        try:
            current_char = nwk_str[char_index]
            if automaton_state == 0:
                # look for the first opening parenthesis
                if current_char == "(":
                    add_new_tree_level()
                    automaton_state = 1
            elif automaton_state == 1 or automaton_state == 3:
                # case 1: name
                # case 3: branch length
                if current_char == ":":
                    automaton_state = 3
                elif current_char == "," or current_char == ")":
                    try:
                        finish_node_definition()
                        automaton_state = 1
                        if current_char == ",":
                            add_new_tree_level()
                    except Exception:
                        return generate_error(char_index)

                elif current_char == "(":
                    if len(current_node_name) > 0:
                        return generate_error(char_index)
                    else:
                        add_new_tree_level()

                elif current_char in name_quotes:
                    if (
                        automaton_state == 1
                        and len(current_node_name) == 0
                        and len(current_node_attribute) == 0
                        and len(current_node_annotation) == 0
                    ):
                        automaton_state = 2
                        quote_delimiter = current_char
                        continue
                    return generate_error(char_index)
                else:
                    if current_char == "{":
                        if len(current_node_annotation):
                            return generate_error(char_index)
                        else:
                            automaton_state = 4
                    else:
                        if automaton_state == 3:
                            current_node_attribute += current_char
                        else:
                            if space.search(current_char):
                                continue
                            if current_char == ";":
                                char_index = len(nwk_str)
                                break
                        current_node_name += current_char
            elif automaton_state == 2:
                # inside a quoted expression
                if current_char == quote_delimiter:
                    if char_index < len(nwk_str) - 1:
                        if nwk_str[char_index + 1] == quote_delimiter:
                            char_index += 1
                            current_node_name += quote_delimiter
                            continue

                    quote_delimiter = 0
                    automaton_state = 1
                    continue
                else:
                    current_node_name += current_char
            elif automaton_state == 4:
                # inside a comment / attribute
                if current_char == "}":
                    automaton_state = 3
                else:
                    if current_char == "{":
                        return generate_error(char_index)
                    current_node_annotation += current_char
        except Exception:
            return generate_error(char_index)

    if len(clade_stack) != 1:
        return generate_error(len(nwk_str) - 1)

    if len(current_node_name):
        if bootstrap_values and "children" in tree_json:
            tree_json["bootstrap_values"] = current_node_name
        else:
            tree_json["name"] = current_node_name

    return tree_json


def traverse_tree(
    node: Dict,
    parent: Optional[Dict],
    labels: Dict,
    labeler: Dict,
    composition: Dict,
    subs: Dict,
    leaf_label: Optional[str] = None,
    ignore_leaves: bool = False,
) -> None:
    """Traverse a phylogenetic tree and collect information about compositions and substitutions.

    Args:
        node: Current node in the tree
        parent: Parent node (None for root)
        labels: Node labels with codon information
        labeler: Dictionary mapping node names to tags
        composition: Dictionary to store composition information
        subs: Dictionary to store substitution information
        leaf_label: Label to use for leaf nodes. Takes precedence over ignore_leaves.
        ignore_leaves: Boolean, whether to calculate composition from leaf node parents' tags (not leaf tags). Useful when tagging internal branches only.
    """
    tag = labeler[node["name"]] if parent else None

    # Track parent-child relationships for substitution tracking
    if parent and "tag" in parent and "tag" in node and parent["tag"] != node["tag"]:
        transition_key = f"{parent['tag']}->{node['tag']}"
        if transition_key not in subs:
            subs[transition_key] = Counter()
        # Add a dummy counter to ensure the transition is tracked
        subs[transition_key]["transition"] = 1

    if node["name"] in labels:
        node["label"] = labels[
            node["name"]
        ]  # what codon is present at the current node?
        if parent and parent["label"]:
            diff = 0
            for i, c in enumerate(node["label"]):
                if c in nucs and i < len(parent["label"]):
                    if parent["label"][i] != c and parent["label"][i] in nucs:
                        diff += 1  # count number of substitutions between parent node and child node

            if diff > 0:
                if tag not in subs:
                    subs[tag] = Counter()

                pt = TT(parent["label"])
                nt = TT(node["label"])
                if pt < nt:
                    sub = pt + ":" + nt
                else:
                    sub = nt + ":" + pt

                subs[tag][sub] += 1

    else:
        node["label"] = parent["label"] if parent else ""

    if "children" not in node:
        if leaf_label:
            tag = leaf_label
        elif ignore_leaves:
            tag = parent["tag"] if parent and "tag" in parent else tag
        else:
            tag = tag

        if tag not in composition:
            composition[tag] = Counter()
        composition[tag][TT(node["label"])] += 1

    if "children" in node:
        for c in node["children"]:
            traverse_tree(
                c, node, labels, labeler, composition, subs, leaf_label, ignore_leaves
            )
