"""
Tree manipulation utilities for HyPhy results.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

def traverse_tree(node: dict, parent: dict, labels: dict, labeler: callable, composition: dict, subs: dict) -> None:
    """Traverse a phylogenetic tree and collect information about compositions and substitutions.
    
    Args:
        node: Current node in the tree
        parent: Parent node
        labels: Node labels
        labeler: Function to determine node label
        composition: Dictionary to store composition information
        subs: Dictionary to store substitution information
    """
    if parent is not None:
        my_label = labeler(node)
        if my_label not in composition:
            composition[my_label] = {}
            subs[my_label] = {}
            
        if node["name"] in labels:
            for k, v in labels[node["name"]].items():
                if k not in composition[my_label]:
                    composition[my_label][k] = 0
                composition[my_label][k] += v
                
                if k not in subs[my_label]:
                    subs[my_label][k] = 0
                subs[my_label][k] += 1
                
    if "children" in node:
        for child in node["children"]:
            traverse_tree(child, node, labels, labeler, composition, subs)
