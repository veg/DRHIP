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
    'prime_marker',
]

COMPARISON_GROUPS_SUMMARY_FIELDNAMES = [
    'gene',
    'comparison_group',
    'group_N',
    'group_T',
    'group_dN/dS',
    'group_nt_conserved',
    'group_aa_conserved',
]

COMPARISON_GROUPS_SITE_FIELDNAMES = [
    'gene',
    'site',
    'comparison_group',
    'unique_aa',
    'diff_majority_residue',
    'cfel_marker',
    'intensified_positive_selection'
]
