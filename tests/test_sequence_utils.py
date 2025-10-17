"""
Tests for sequence utilities.
"""

from collections import Counter

import pytest

from drhip.utils import sequence_utils as su


def test_get_majority_residue():
    """Test getting the majority residue from a composition."""
    # Test with normal composition
    composition = {"A": 5, "G": 3, "T": 1}
    assert su.get_majority_residue(composition) == "A"

    # Test with empty composition
    assert su.get_majority_residue({}) == "-"

    # Test with tied frequencies
    composition = {"A": 3, "G": 3, "T": 1}
    # Should return either A or G (depends on implementation)
    result = su.get_majority_residue(composition)
    assert result in ["A", "G"]


def test_has_diff_majority_residue():
    """Test checking if majority residue differs between clades."""
    # Test with different majority residue in second clade
    focal = Counter({"A": 5, "G": 2})
    others = [Counter({"A": 4, "G": 3}), Counter({"G": 5, "A": 2})]
    assert su.has_diff_majority_residue(focal, others)
    # Test with same majority residue
    focal = Counter({"A": 5, "G": 2})
    others = [Counter({"A": 4, "G": 3})]
    assert not su.has_diff_majority_residue(focal, others)

    # Test with empty compositions
    assert not su.has_diff_majority_residue({}, [])
    assert not su.has_diff_majority_residue(focal, [])
    assert not su.has_diff_majority_residue({}, others)

    # Test with multiple other clades
    focal = Counter({"A": 5, "G": 2})
    others = [
        Counter({"A": 4, "G": 3}),  # Same majority
        Counter({"G": 4, "A": 1}),  # Different majority
    ]
    assert su.has_diff_majority_residue(focal, others)


def test_get_unique_aa():
    """Test getting amino acids unique to the focal clade."""
    # Test with unique amino acids
    focal = Counter({"A": 5, "G": 2, "T": 1})
    others = [Counter({"A": 4, "G": 3})]
    result = su.get_unique_aa(focal, others)
    assert "T" in result
    assert "A" not in result
    assert "G" not in result

    # Test with no unique amino acids
    focal = Counter({"A": 5, "G": 2})
    others = [Counter({"A": 4, "G": 3, "T": 1})]
    assert su.get_unique_aa(focal, others) == ""

    # Test with empty compositions
    assert su.get_unique_aa({}, []) == ""
    assert su.get_unique_aa(focal, []) == "A G"  # All are unique if no other clades

    # Test with multiple unique amino acids
    focal = Counter({"A": 5, "G": 2, "T": 1, "C": 1})
    others = [Counter({"A": 4})]
    result = su.get_unique_aa(focal, others)
    assert set(result.split()) == {"G", "T", "C"}


def test_format_composition():
    """Test formatting composition data."""
    # Test with normal composition
    composition = Counter({"A": 5, "G": 3, "T": 1})
    result = su.format_composition(composition)
    assert "A:5" in result
    assert "G:3" in result
    assert "T:1" in result

    # Test with empty composition
    assert su.format_composition({}) == "-"


def test_format_substitutions():
    """Test formatting substitution data."""
    # Test with normal substitutions
    subs = Counter({"A->G": 5, "T->C": 2})
    result = su.format_substitutions(subs)
    assert "A->G:5" in result
    assert "T->C:2" in result

    # Test with empty substitutions
    assert su.format_substitutions({}) == "-"


def test_get_site_composition():
    """Test calculating site composition from sequences."""
    # Test with normal sequences
    sequences = {"seq1": "AGCT", "seq2": "AGTT", "seq3": "AGGT"}

    # Test site 0 (all A)
    result = su.get_site_composition(sequences, 0)
    assert result == {"A": 1.0}

    # Test site 1 (all G)
    result = su.get_site_composition(sequences, 1)
    assert result == {"G": 1.0}

    # Test site 2 (mixed C and G)
    result = su.get_site_composition(sequences, 2)
    assert abs(result["C"] - 1 / 3) < 0.001
    assert abs(result["G"] - 1 / 3) < 0.001

    # Test site 3 (all T)
    result = su.get_site_composition(sequences, 3)
    assert result == {"T": 1.0}


def test_get_majority_residue_from_frequencies():
    """Test getting majority residue from frequency dictionary."""
    # Test with normal frequencies
    composition = {"A": 0.6, "G": 0.3, "T": 0.1}
    assert su.get_majority_residue_from_frequencies(composition) == "A"

    # Test with empty composition
    assert su.get_majority_residue_from_frequencies({}) == "-"

    # Test with tied frequencies
    composition = {"A": 0.4, "G": 0.4, "T": 0.2}
    # Should return either A or G (depends on implementation)
    result = su.get_majority_residue_from_frequencies(composition)
    assert result in ["A", "G"]


def test_get_unique_aa_count():
    """Test counting unique amino acids."""
    # Test with normal composition
    composition = {"A": 0.5, "G": 0.3, "T": 0.2}
    assert su.get_unique_aa_count(composition) == 3

    # Test with zero frequencies
    composition = {"A": 0.8, "G": 0.2, "T": 0.0}
    assert su.get_unique_aa_count(composition) == 2

    # Test with empty composition
    assert su.get_unique_aa_count({}) == 0


def test_compare_majority_residues():
    """Test comparing majority residues between compositions."""
    # Test with different majority residues
    reference = {"A": 0.6, "G": 0.4}
    target = {"G": 0.7, "A": 0.3}
    assert su.compare_majority_residues(reference, target)

    # Test with same majority residue
    reference = {"A": 0.6, "G": 0.4}
    target = {"A": 0.8, "G": 0.2}
    assert not su.compare_majority_residues(reference, target)

    # Test with empty compositions
    assert not su.compare_majority_residues({}, {})
    assert not su.compare_majority_residues(reference, {})
    assert not su.compare_majority_residues({}, target)


def test_process_sequence_data():
    """Test processing sequence data for site-specific metrics."""
    # Create test sequence data
    sequences_by_group = {
        "group1": {"seq1": "AGCT", "seq2": "AGCT"},
        "group2": {"seq3": "AGTT", "seq4": "AGGT"},
    }

    # Process all sites
    sites = [0, 1, 2, 3]
    result = su.process_sequence_data(sequences_by_group, sites)

    # Check that we have data for all sites
    assert set(result.keys()) == set(sites)

    # Check site 0 (all A in both groups)
    assert not result[0]["diff_majority_residue"]

    # Check site 2 (C in group1, G in group2)
    assert result[2]["diff_majority_residue"]

    # Check that each site has the expected metrics
    for site_idx, site_data in result.items():
        assert "diff_majority_residue" in site_data
        assert "unique_aa" in site_data
        assert "group1_composition" in site_data
        assert "group2_composition" in site_data
