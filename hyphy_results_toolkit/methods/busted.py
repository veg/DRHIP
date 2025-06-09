"""
BUSTED (Branch-Site Unrestricted Statistical Test for Episodic Diversification) method implementation.
"""

from typing import Dict, Any, List

from ..utils import tree_helpers
from .base import HyPhyMethod

class BustedMethod(HyPhyMethod):
    """Implementation of BUSTED analysis processing."""
    
    def __init__(self):
        """Initialize BUSTED method."""
        super().__init__("BUSTED", "BUSTED.json")
    
    def get_omega3(self, fit: Dict[str, Any]) -> Dict[str, float]:
        """Get the omega3 value and proportion from a model fit.
        
        Args:
            fit: Dictionary containing model fit results
            
        Returns:
            Dictionary containing omega and proportion values
        """
        omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
        omega_idx = str(len(omegas) - 1)
        return omegas[omega_idx]
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process BUSTED results.
        
        Args:
            results: Raw BUSTED results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        # Extract common fields (N, T, sites)
        processed = self.extract_common_fields(results)
        
        # Add BUSTED-specific fields
        test_results = results['test results']
        processed['BUSTED_pval'] = test_results['p-value']  # Only keep the p-value field
        
        # Get omega3 distribution
        omega3 = self.get_omega3(results)
        processed['BUSTED_omega3'] = omega3['omega']  # Renamed to match desired format
        
        # Get proportion of sites in omega3 category and convert to percentage
        if 'weight' in omega3:
            processed['BUSTED_prop_sites_in_omega3'] = omega3['weight'] * 100.0  # Convert to percentage
        elif 'proportion' in omega3:
            processed['BUSTED_prop_sites_in_omega3'] = omega3['proportion'] * 100.0  # Convert to percentage
        else:
            processed['BUSTED_prop_sites_in_omega3'] = 0.0
        
        # Calculate dN/dS using the helper function
        dnds_stats = self.calculate_rate_distribution_stats(
            results, 
            model_name='Unconstrained model', 
            distribution_name='global'
        )
        
        # Update processed results with dN/dS
        processed.update(dnds_stats)
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
        """Process site-specific data from BUSTED results.
        
        Args:
            results: Raw BUSTED results dictionary
            
        Returns:
            Dictionary mapping site indices to site-specific data
        """
        site_data = {}
        
        # Check if we have substitutions data
        if 'substitutions' not in results or '0' not in results['substitutions']:
            return site_data
            
        # Get the substitutions data
        substitutions = results['substitutions']['0']
        
        # Get the tree
        if 'input' not in results or 'trees' not in results['input'] or '0' not in results['input']['trees']:
            return site_data
            
        # Parse the tree
        internal_branches = {}
        tree = tree_helpers.newick_parser(results['input']['trees']['0'], {}, internal_branches)
        
        if tree['error'] is not None:
            print(f"Error parsing tree: {tree['error']}")
            return site_data
            
        tree = tree['json']
        
        # Process each site
        for site_idx, site_subs in substitutions.items():
            site_num = int(site_idx)
            
            # Initialize composition and substitution counters
            composition = {}
            subs = {}
            
            # Check if we have comparison groups available
            has_comparison_groups = False
            comparison_groups = []
            if hasattr(self, '_comparison_groups') and self._comparison_groups:
                if len(self._comparison_groups) > 0:
                    has_comparison_groups = True
                    comparison_groups = self._comparison_groups
            
            # Initialize composition and substitution dictionaries
            for group in comparison_groups:
                # Traverse the tree to collect composition and substitution data for each group
                tree_helpers.traverse_tree(
                    tree, 
                    None, 
                    site_subs, 
                    internal_branches, 
                    composition, 
                    subs, 
                    group  # Use each group as leaf label
                )
            
            # Process the composition data
            site_composition = []
            for tag, counts in composition.items():
                for aa, count in counts.items():
                    site_composition.append(f"{aa}:{count}")
            
            # Process the substitution data
            site_substitutions = []
            for tag, counts in subs.items():
                for sub, count in counts.items():
                    site_substitutions.append(f"{sub}:{count}")
            
            # Initialize site data with basic information - only include site-specific fields
            site_info = {
                'composition': ','.join(site_composition) if site_composition else 'NA',
                'substitutions': ','.join(site_substitutions) if site_substitutions else 'NA',
                'majority_residue': 'NA'
            }
            
            # Calculate majority residue if we have comparison groups
            if has_comparison_groups and comparison_groups[0] in composition and composition[comparison_groups[0]]:
                sorted_residues = sorted(
                    [[aa, count] for aa, count in composition[comparison_groups[0]].items()],
                    key=lambda d: -d[1]
                )
                if sorted_residues:
                    site_info['majority_residue'] = sorted_residues[0][0]
            
            # Store the site data
            site_data[site_num] = site_info
            
        return site_data
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'BUSTED_pval',
            'BUSTED_omega3',
            'BUSTED_prop_sites_in_omega3'
        ]
        
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'composition',
            'substitutions',
            'majority_residue'
        ]
        
    def get_comparison_group_fields(self, results: Dict[str, Any] = None) -> List[str]:
        """Get list of fields that are specific to comparison groups.
        
        Args:
            results: Optional raw BUSTED results dictionary to check if data is available
            
        Returns:
            List of comparison group fields if requirements are met, empty list otherwise
        """
        # Basic check: return empty list if no comparison groups are set
        if not self._comparison_groups:
            return []
            
        # If results are provided, perform thorough checks
        if results is not None:
            # Check if we have MLE content
            if not self.has_mle_content(results):
                return []
                
            # Check if we have substitutions data
            if 'substitutions' not in results or '0' not in results['substitutions']:
                return []
                
            # Check if we have tree data
            if 'input' not in results or 'trees' not in results['input'] or '0' not in results['input']['trees']:
                return []
        
        # If we pass all checks or no results provided, return the fields
        return [
            'diff_majority_residue',
            'unique_aa'
        ]
        
    def process_comparison_site_data(self, results: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from BUSTED results.
        
        Args:
            results: Raw BUSTED results dictionary
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # Skip processing if no comparison groups or no MLE content
        if not self._comparison_groups or not self.has_mle_content(results):
            return {}
            
        # Initialize result dictionary
        comparison_data = {}
        
        # Check if we have substitutions data
        if 'substitutions' not in results or '0' not in results['substitutions']:
            return comparison_data
            
        # Get the substitutions data
        substitutions = results['substitutions']['0']
        
        # Get the tree
        if 'input' not in results or 'trees' not in results['input'] or '0' not in results['input']['trees']:
            return comparison_data
            
        # Parse the tree
        internal_branches = {}
        tree = tree_helpers.newick_parser(results['input']['trees']['0'], {}, internal_branches)
        
        if tree['error'] is not None:
            print(f"Error parsing tree: {tree['error']}")
            return comparison_data
            
        tree = tree['json']
        
        # Process each site
        for site_idx, site_subs in substitutions.items():
            site_num = int(site_idx)
            
            # Initialize composition and substitution counters
            composition = {}
            subs = {}
            
            # Initialize composition and substitution dictionaries
            for group in self._comparison_groups:
                # Traverse the tree to collect composition and substitution data for each group
                tree_helpers.traverse_tree(
                    tree, 
                    None, 
                    site_subs, 
                    internal_branches, 
                    composition, 
                    subs, 
                    group  # Use each group as leaf label
                )
            
            # Initialize site comparison data
            site_comparison_data = {}
            
            # Initialize unique amino acids and different majority residue
            unique_aas_by_group = {}
            diff_majority_residue = False
            
            # For each group, compare against all others
            for i, focal_group in enumerate(self._comparison_groups):
                if focal_group in composition and composition[focal_group]:
                    # Get all amino acids from focal group
                    focal_aas = set(composition[focal_group].keys())
                    
                    # Get majority residue in focal group
                    focal_sorted = sorted(
                        [[aa, count] for aa, count in composition[focal_group].items()],
                        key=lambda d: -d[1]
                    )
                    
                    if focal_sorted:
                        focal_majority = focal_sorted[0][0]
                        
                        # Get all amino acids from other groups
                        other_aas = set()
                        for j, other_group in enumerate(self._comparison_groups):
                            if i != j and other_group in composition:
                                other_aas.update(composition[other_group].keys())
                                
                                # Check for different majority residue
                                if composition[other_group]:
                                    other_sorted = sorted(
                                        [[aa, count] for aa, count in composition[other_group].items()],
                                        key=lambda d: -d[1]
                                    )
                                    
                                    if other_sorted and focal_majority != other_sorted[0][0]:
                                        diff_majority_residue = True
                        
                        # Find unique amino acids for this group
                        unique_aas = focal_aas - other_aas
                        if unique_aas:
                            unique_aas_by_group[focal_group] = ' '.join(sorted(unique_aas))
            
            # Format unique amino acids as a string
            unique_aa_str = ''
            if unique_aas_by_group:
                unique_aa_parts = []
                for group, aas in unique_aas_by_group.items():
                    unique_aa_parts.append(f"{group}:{aas}")
                unique_aa_str = ','.join(unique_aa_parts)
            
            # Add comparison group-specific data for each group
            for group in self._comparison_groups:
                group_data = {
                    'diff_majority_residue': diff_majority_residue,
                    'unique_aa': unique_aa_str if unique_aa_str else 'NA'
                }
                site_comparison_data[group] = group_data
            
            # Only add sites with comparison data
            if site_comparison_data:
                comparison_data[site_num] = site_comparison_data
        
        return comparison_data
