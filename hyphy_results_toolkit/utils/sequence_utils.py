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
    """Get the most frequent amino acid at a site.
    
    Args:
        composition: Dictionary of amino acid counts at a site
        
    Returns:
        The most frequent amino acid or '-' if no data
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
    
    Args:
        focal_composition: Counter of amino acid frequencies in the focal clade
        other_compositions: List of Counters for other clades
        
    Returns:
        True if the majority residue differs, False otherwise
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
    Calculate the amino acid composition at a specific site.
    
    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        site_index: Zero-based index of the site to analyze
        
    Returns:
        Dictionary mapping amino acids to their frequency
    """
    # Extract the residue at the specified site from each sequence
    residues = [seq[site_index] for seq in sequences.values() if site_index < len(seq)]
    
    # Count occurrences of each residue
    residue_counts = Counter(residues)
    
    # Calculate frequencies
    total = sum(residue_counts.values())
    frequencies = {residue: count / total for residue, count in residue_counts.items()}
    
    return frequencies


def get_majority_residue_from_frequencies(composition: Dict[str, float]) -> str:
    """
    Determine the majority residue from a composition dictionary with frequencies.
    
    Args:
        composition: Dictionary mapping residues to their frequencies
        
    Returns:
        The most common residue
    """
    if not composition:
        return "-"
    
    # Find the residue with the highest frequency
    return max(composition.items(), key=lambda x: x[1])[0]


def get_unique_aa_count(composition: Dict[str, float]) -> int:
    """
    Count the number of unique amino acids at a site.
    
    Args:
        composition: Dictionary mapping residues to their frequencies
        
    Returns:
        Number of unique amino acids
    """
    # Count residues with non-zero frequency
    return len([aa for aa, freq in composition.items() if freq > 0])


def compare_majority_residues(reference_composition: Dict[str, float], 
                             target_composition: Dict[str, float]) -> bool:
    """
    Compare majority residues between reference and target compositions.
    
    Args:
        reference_composition: Dictionary mapping residues to frequencies in reference group
        target_composition: Dictionary mapping residues to frequencies in target group
        
    Returns:
        True if majority residues are different, False otherwise
    """
    ref_majority = get_majority_residue(reference_composition)
    target_majority = get_majority_residue(target_composition)
    
    return ref_majority != target_majority


def process_sequence_data(sequences_by_group: Dict[str, Dict[str, str]], 
                         sites: List[int]) -> Dict[int, Dict[str, any]]:
    """
    Process sequence data to extract site-specific metrics.
    
    Args:
        sequences_by_group: Dictionary mapping group names to sequence dictionaries
        sites: List of site indices to process
        
    Returns:
        Dictionary mapping site indices to dictionaries of site-specific metrics
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
