"""
Utilities for sequence analysis and site-specific data processing.

This module provides functions for analyzing sequence data and generating
site-specific metrics that are not directly provided by HyPhy methods.

Sequence analysis utilities for HyPhy results.

Authors:
    Danielle Callan (dcallan@temple.edu)
"""

from typing import Dict, List
from collections import Counter


def get_majority_residue(composition: Dict[str, Counter]) -> str:
    """Get the most frequent amino acid at a site based on count data.
    
    This function identifies the amino acid with the highest count in the provided
    composition dictionary. In case of ties, it returns the first one based on
    sorting order.
    
    Args:
        composition: Dictionary mapping amino acids to their counts at a specific site
        
    Returns:
        The most frequent amino acid or '-' if the composition is empty
    """
    if not composition:
        return '-'
    
    # Sort by frequency (descending) and return the most common residue
    sorted_residues = sorted(
        [[aa, count] for aa, count in composition.items()], 
        key=lambda d: -d[1]
    )
    
    return sorted_residues[0][0] if sorted_residues else '-'


def has_diff_majority_residue(focal_composition: Counter, other_compositions: List[Counter]) -> bool:
    """Check if the majority residue differs between the focal clade and any other clade.
    
    This function compares the most common amino acid in the focal composition with
    the most common amino acid in each of the other compositions. If any of them differ,
    it returns True. This is useful for identifying sites where evolutionary pressures
    may have caused different amino acids to dominate in different clades.
    
    Args:
        focal_composition: Counter of amino acid frequencies in the focal clade
        other_compositions: List of Counters for other clades to compare against
        
    Returns:
        True if the majority residue in the focal clade differs from any other clade,
        False if they all match or if either input is empty
    """
    if not focal_composition or not other_compositions:
        return False
    
    # Get the majority residue in the focal clade
    focal_cons = sorted(
        [[aa, count] for aa, count in focal_composition.items()], 
        key=lambda d: -d[1]
    )
    
    if not focal_cons:
        return False
    
    focal_majority = focal_cons[0][0]
    
    # Compare with each other clade
    for other_comp in other_compositions:
        if not other_comp:
            continue
            
        other_cons = sorted(
            [[aa, count] for aa, count in other_comp.items()], 
            key=lambda d: -d[1]
        )
        
        if not other_cons:
            continue
            
        if focal_majority != other_cons[0][0]:
            return True
    
    return False


def get_unique_aa(focal_composition: Counter, other_compositions: List[Counter]) -> str:
    """Get amino acids unique to the focal clade.
    
    Args:
        focal_composition: Counter of amino acid frequencies in the focal clade
        other_compositions: List of Counters for other clades
        
    Returns:
        Space-separated string of unique amino acids or empty string if none
    """
    if not focal_composition:
        return ''
    
    # Collect all amino acids from other clades
    other_clade_aas = set()
    for comp in other_compositions:
        other_clade_aas.update(comp.keys())
    
    # Find amino acids unique to the focal clade
    unique_aas = set(focal_composition.keys()) - other_clade_aas
    
    return ' '.join(sorted(unique_aas))


def format_composition(composition: Counter) -> str:
    """Format composition data for output.
    
    Args:
        composition: Counter of amino acid frequencies
        
    Returns:
        Formatted string representation of composition
    """
    if not composition:
        return '-'
    
    return ', '.join([f'{aa}:{count}' for aa, count in composition.items()])


def format_substitutions(substitutions: Counter) -> str:
    """Format substitution data for output.
    
    Args:
        substitutions: Counter of substitution types
        
    Returns:
        Formatted string representation of substitutions
    """
    if not substitutions:
        return '-'
    
    return ', '.join([f'{sub}:{count}' for sub, count in substitutions.items()])


def get_site_composition(sequences: Dict[str, str], site_index: int) -> Dict[str, float]:
    """
    Calculate the amino acid composition at a specific site across multiple sequences.
    
    This function extracts the residue at the specified site from each sequence,
    counts the occurrences of each residue, and calculates their frequencies.
    It's a fundamental function for analyzing site-specific evolutionary patterns.
    
    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        site_index: Zero-based index of the site to analyze
        
    Returns:
        Dictionary mapping amino acids to their frequency (values sum to 1.0)
    """
    # Extract the residue at the specified site from each sequence
    residues = [seq[site_index] for seq in sequences.values() if site_index < len(seq)]
    
    # Count occurrences of each residue
    residue_counts = Counter(residues)
    
    # Calculate frequencies
    total = sum(residue_counts.values())
    frequencies = {residue: count / total for residue, count in residue_counts.items()}
    
    return frequencies


def get_majority_residue_from_frequencies(frequencies: Dict[str, float]) -> str:
    """
    Get the majority residue from a frequency dictionary.
    
    Similar to get_majority_residue but works with frequency data (0.0-1.0) rather than
    raw counts. This function is particularly useful when working with normalized data
    or when comparing compositions across different sample sizes.
    
    Args:
        frequencies: Dictionary mapping residues to their frequencies (values between 0.0 and 1.0)
        
    Returns:
        The residue with the highest frequency or '-' if the dictionary is empty
    """
    if not frequencies:
        return "-"
    
    # Find the residue with the highest frequency
    return max(frequencies.items(), key=lambda x: x[1])[0]


def get_unique_aa_count(composition: Dict[str, float]) -> int:
    """
    Count the number of unique amino acids in a composition.
    
    This function simply returns the number of different amino acids present at a site,
    regardless of their frequencies. It's a useful metric for measuring site diversity
    and can indicate sites under different types of selection pressure.
    
    Args:
        composition: Dictionary mapping amino acids to their frequencies
        
    Returns:
        Number of unique amino acids (integer count of dictionary keys)
    """
    # Count residues with non-zero frequency
    return len([aa for aa, freq in composition.items() if freq > 0])


def compare_majority_residues(reference_composition: Dict[str, float], 
                             target_composition: Dict[str, float]) -> bool:
    """
    Compare majority residues between reference and target compositions.
    
    Unlike has_diff_majority_residue which works with multiple compositions,
    this function performs a direct one-to-one comparison between two frequency
    dictionaries. It uses get_majority_residue_from_frequencies to extract the
    majority residue from each composition.
    
    Args:
        reference_composition: Dictionary mapping residues to frequencies in reference group
        target_composition: Dictionary mapping residues to frequencies in target group
        
    Returns:
        True if majority residues are different, False if they match or if either input is empty
    """
    # If either composition is empty, return False
    if not reference_composition or not target_composition:
        return False
        
    ref_majority = get_majority_residue_from_frequencies(reference_composition)
    target_majority = get_majority_residue_from_frequencies(target_composition)
    
    return ref_majority != target_majority


def process_sequence_data(sequences_by_group: Dict[str, Dict[str, str]], 
                         sites: List[int]) -> Dict[int, Dict[str, any]]:
    """
    Process sequence data to extract site-specific metrics across different groups.
    
    This function analyzes sequence data organized by groups (e.g., different viral strains
    or evolutionary clades) and calculates various metrics for each specified site, including:
    - Composition of amino acids for each group
    - Overall majority residue across all sequences
    - Count of unique amino acids at the site
    - Whether the majority residue differs between the reference group and other groups
    
    The first group in sequences_by_group is considered the reference group for comparisons.
    
    Args:
        sequences_by_group: Dictionary mapping group names to dictionaries of sequence IDs and sequences
        sites: List of site indices to process (zero-based)
        
    Returns:
        Dictionary mapping site indices to dictionaries containing site-specific metrics
        including group-specific compositions and comparison results
    """
    result = {}
    
    # Get reference group if available (usually the first group)
    reference_group = next(iter(sequences_by_group)) if sequences_by_group else None
    
    for site in sites:
        site_data = {}
        
        # Get composition for each group
        compositions = {}
        for group, sequences in sequences_by_group.items():
            compositions[group] = get_site_composition(sequences, site)
            # Add group-specific composition to site data
            site_data[f'{group}_composition'] = compositions[group]
        
        # Calculate overall composition (combine all groups)
        all_sequences = {}
        for group_sequences in sequences_by_group.values():
            all_sequences.update(group_sequences)
        
        overall_composition = get_site_composition(all_sequences, site)
        
        # Get majority residue
        site_data['majority_residue'] = get_majority_residue(overall_composition)
        
        # Get unique AA count
        site_data['unique_aa'] = get_unique_aa_count(overall_composition)
        
        # Check if majority residue differs between reference and other groups
        if reference_group and len(sequences_by_group) > 1:
            other_groups = {k: v for k, v in sequences_by_group.items() if k != reference_group}
            
            # Combine all other groups
            other_sequences = {}
            for group_sequences in other_groups.values():
                other_sequences.update(group_sequences)
            
            other_composition = get_site_composition(other_sequences, site)
            
            # Check if majority residue differs
            site_data['diff_majority_residue'] = compare_majority_residues(
                compositions[reference_group], other_composition
            )
        else:
            # If there's only one group, there can't be a difference
            site_data['diff_majority_residue'] = False
        
        result[site] = site_data
    
    return result
