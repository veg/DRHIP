"""
Tests for tree helper functions.
"""

import pytest
from collections import Counter
from unittest.mock import patch, MagicMock

from drhip.utils import tree_helpers as th


def test_tt_translation():
    """Test codon translation function."""
    # Test valid codons
    assert th.TT("ATG") == "M"  # Methionine
    assert th.TT("TAA") == "*"  # Stop codon
    assert th.TT("GGG") == "G"  # Glycine
    
    # Test invalid codon
    assert th.TT("XYZ") == "?"
    
    # Test gap codon
    assert th.TT("---") == "-"


def test_newick_parser_basic():
    """Test basic Newick tree parsing."""
    # Simple tree
    newick_str = "(A,B);"
    result = th.newick_parser(newick_str, False)
    
    # Check structure
    assert "json" not in result  # No error
    assert "name" in result
    assert "children" in result
    assert len(result["children"]) == 2
    
    # Check node names
    child_names = [child["name"] for child in result["children"]]
    assert "A" in child_names
    assert "B" in child_names


def test_newick_parser_with_bootstrap():
    """Test Newick parsing with bootstrap values."""
    # Tree with bootstrap values
    newick_str = "(A,B)0.95;"
    result = th.newick_parser(newick_str, True)
    
    # Check structure
    assert "json" not in result  # No error
    assert "bootstrap_values" in result
    assert result["bootstrap_values"] == "0.95"


def test_newick_parser_with_branch_lengths():
    """Test Newick parsing with branch lengths."""
    # Tree with branch lengths
    newick_str = "(A:0.1,B:0.2);"
    result = th.newick_parser(newick_str, False)
    
    # Check structure
    assert "json" not in result  # No error
    
    # Check branch lengths (stored in attribute)
    child_attributes = [child["attribute"] for child in result["children"]]
    assert "0.1" in child_attributes
    assert "0.2" in child_attributes


def test_newick_parser_with_quotes():
    """Test Newick parsing with quoted node names."""
    # Tree with quoted names
    newick_str = "('Node A','Node B');"
    result = th.newick_parser(newick_str, False)
    
    # Check node names
    child_names = [child["name"] for child in result["children"]]
    assert "Node A" in child_names
    assert "Node B" in child_names


def test_newick_parser_with_tags():
    """Test Newick parsing with node tags."""
    # Tree with tags to track
    newick_str = "(A,B);"
    track_tags = {}
    optional_tags = {"A": "foreground", "B": "background"}
    
    result = th.newick_parser(newick_str, False, track_tags, optional_tags)
    
    # Check tags were tracked
    assert "A" in track_tags
    assert track_tags["A"] == "foreground"
    assert "B" in track_tags
    assert track_tags["B"] == "background"


def test_newick_parser_nested():
    """Test Newick parsing with nested structure."""
    # Nested tree
    newick_str = "((A,B),(C,D));"
    result = th.newick_parser(newick_str, False)
    
    # Check structure
    assert "json" not in result  # No error
    assert len(result["children"]) == 2
    
    # Check that both children have children
    for child in result["children"]:
        assert "children" in child
        assert len(child["children"]) == 2


def test_newick_parser_error():
    """Test Newick parsing with invalid input."""
    # Invalid Newick string (unbalanced parentheses)
    newick_str = "((A,B),(C,D);"
    result = th.newick_parser(newick_str, False)
    
    # Check error
    assert "json" in result
    assert result["json"] is None
    assert "error" in result


def test_traverse_tree():
    """Test tree traversal function."""
    # Create a simple tree
    tree = {
        "name": "root",
        "children": [
            {
                "name": "node1",
                "tag": "foreground"
            },
            {
                "name": "node2",
                "tag": "background",
                "children": [
                    {
                        "name": "node3",
                        "tag": "background"
                    }
                ]
            }
        ]
    }
    
    # Create mock data
    labels = {
        "node1": "ATG",
        "node2": "ATG",
        "node3": "ATT"
    }
    
    labeler = {
        "node1": "foreground",
        "node2": "background",
        "node3": "background"
    }
    
    # Initialize tracking dictionaries
    composition = {}
    subs = {}
    
    # Traverse tree
    th.traverse_tree(tree, None, labels, labeler, composition, subs)
    
    # Check composition tracking
    assert "foreground" in composition
    assert "background" in composition
    assert composition["foreground"] == {"M": 1}
    assert composition["background"] == {"I": 1}
    
    # Check substitution tracking - should have background tag with substitutions
    assert "background" in subs
    assert "I:M" in subs["background"] or "M:I" in subs["background"]

def test_traverse_tree_group_specific():
    """Test tree traversal function."""
    # Create a simple tree
    tree = {
        "name": "root",
        "children": [
            {
                "name": "node1",
                "tag": "foreground",
                "children": [
                    {
                        "name": "node2",
                        "tag": "background"
                    }
                ]
            },
            {
                "name": "node3",
                "tag": "reference",
                "children": [
                    {
                        "name": "node4",
                        "tag": "background"
                    }
                ]
            }
        ]
    }
    
    # Create mock data
    labels = {
        "node1": "ATG",
        "node2": "ATG",
        "node3": "ATT",
        "node4": "ATT"
    }
    
    labeler = {
        "node1": "foreground",
        "node2": "background",
        "node3": "reference",
        "node4": "background"
    }
    
    # Initialize tracking dictionaries
    composition = {}  # fucker only counts leaves
    subs = {}
    
    # Traverse tree
    th.traverse_tree(tree, None, labels, labeler, composition, subs, "foreground")
    
    # Check composition tracking
    assert "foreground" in composition
    assert "background" in composition
    assert "reference" in composition
    assert composition["foreground"] == {"M": 1}
    assert composition["background"] == {"I": 1}
    # assert composition["reference"] == {"M": 1}
    
    # Check substitution tracking - should have background tag with substitutions
    assert "background" in subs
    assert "I:M" in subs["background"] or "M:I" in subs["background"]
