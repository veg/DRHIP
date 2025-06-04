"""
Configuration settings for the hyphy-results-toolkit.
"""

# Default values for comparison groups
DEFAULT_COMPARISON_GROUPS = ['foreground', 'background']

# Method-specific paths relative to results directory for CAPHEINE structure
METHOD_PATHS = {
    'BUSTED': 'BUSTED',
    'RELAX': 'RELAX',
    'CFEL': 'contrastFEL',
    'FEL': 'FEL',
    'MEME': 'MEME',
    'PRIME': 'PRIME'
}

# The values in these lists are populated by individual methods
# The lists here are used to validate that all expected fields are present in the output files
SUMMARY_FIELDNAMES = [
    'gene',
    'N',      # Number of sequences
    'T',      # Total branch length
    'dN/dS',  # Overall dN/dS ratio
    'sites',  # Number of sites
    'diff_sites',
]

SITES_FIELDNAMES = [
    'gene',
    'site',
    'composition',
    'substitutions',
    'majority_residue',
    'diff_majority_residue',
    'unique_aa',
    'intensified_positive_selection',
    'cfel_marker',
    'prime_marker',
]
