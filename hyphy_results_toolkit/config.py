"""
Configuration settings for the hyphy-results-toolkit.
"""

# Method-specific paths relative to results directory for CAPHEINE structure
METHOD_PATHS = {
    'BUSTED': 'BUSTED',
    'RELAX': 'RELAX',
    'CFEL': 'contrastFEL',
    'FEL': 'FEL',
    'MEME': 'MEME',
    'PRIME': 'PRIME'
}

SUMMARY_FIELDNAMES = [
    'gene',
    'comparison_group',  # Keep using comparison_group internally
    'N',      # Number of sequences
    'T',      # Total branch length
    'dN/dS',  # Overall dN/dS ratio
    'sites',  # Number of sites
    'nt_conserved',
    'aa_conserved',
    'diff_sites',
    # Method-specific fields will be added by each method
]

SITES_FIELDNAMES = [
    'gene',
    'comparison_group',  # Keep using comparison_group internally
    'site',
    'consensus_site',
    'composition',
    'substitutions',
    'majority_residue',
    'diff_majority_residue',
    'unique_aa',
    'intensified_positive_selection',
    # Method-specific fields will be added by each method
]
