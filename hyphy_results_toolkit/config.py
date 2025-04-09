"""
Configuration settings for the hyphy-results-toolkit.
"""

# Method-specific paths relative to results directory
METHOD_PATHS = {
    'BUSTED': 'concat/BUSTED',
    'RELAX': 'concat/RELAX',
    'CFEL': 'concat/contrastFEL',
    'FEL': 'concat/FEL',
    'MEME': 'concat/MEME',
    'PRIME': 'concat/PRIME'
}

# Default field names for output files
SUMMARY_FIELDNAMES = [
    'gene',
    'clade',
    'nt_conserved',
    'aa_conserved',
    'total_branch_length',
    'tree_length_ratio',
    # Method-specific fields will be added by each method
]

SITES_FIELDNAMES = [
    'gene',
    'site',
    'clade',
    # Method-specific fields will be added by each method
]
