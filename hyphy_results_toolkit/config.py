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

# Default field names for output files
SUMMARY_FIELDNAMES = [
    'gene',
    'comparison_group',  # Previously 'clade'
    'nt_conserved',
    'aa_conserved',
    'total_branch_length',
    'tree_length_ratio',
    # Method-specific fields will be added by each method
]

SITES_FIELDNAMES = [
    'gene',
    'site',
    'comparison_group',  # Previously 'clade'
    # Method-specific fields will be added by each method
]
