"""
BUSTED (Branch-Site Unrestricted Statistical Test for Episodic Diversification) method implementation.
"""

from typing import Dict, Any, List

from ..utils import tree_helpers
from ..utils import sequence_utils as su
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
        
    def calculate_rate_distribution_stats(self, results: Dict[str, Any],
                                     model_name: str = 'Unconstrained model',
                                     distribution_name: str = 'global') -> Dict[str, float]:
        """Calculate statistics from rate distributions.
        
        This helper extracts rate distributions and calculates statistics like dN/dS.
        
        Args:
            results: Raw results dictionary from JSON file
            model_name: Name of the model to extract rates from (default: 'Unconstrained model')
            distribution_name: Name of the distribution to use (default: 'global')
            
        Returns:
            Dictionary with calculated statistics (e.g., {'dN/dS': 0.5})
        """
        stats = {'dN/dS': 0.0}
        
        try:
            if 'fits' in results and model_name in results['fits']:
                model_fit = results['fits'][model_name]
                if 'Rate Distributions' in model_fit and distribution_name in model_fit['Rate Distributions']:
                    rates = model_fit['Rate Distributions'][distribution_name]
                    
                    # Calculate weighted average of omega values
                    omega_sum = 0.0
                    for rate in rates:
                        # Handle different formats (weight or proportion)
                        weight = rate.get('weight', rate.get('proportion', 0.0))
                        omega_sum += rate['omega'] * weight
                    
                    stats['dN/dS'] = omega_sum
        except Exception as e:
            print(f"Error calculating rate distribution stats: {e}")
            # Keep default values on error
        
        return stats
    
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

        # Process each site
        for site_idx, site_subs in substitutions.items():
            # Convert to 1-based indexing to match FEL and other methods
            site_num = int(site_idx) + 1
            
            # Initialize composition and substitution dictionaries
            composition = {}
            subs = {}
            
            # Traverse the tree to collect composition and substitution data for each group
            tree_helpers.traverse_tree(
                tree, 
                None, 
                site_subs, 
                internal_branches, 
                composition, 
                subs, 
                None  # no leaf labels
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
            
            # Calculate majority residue for test (internal) branches
            if 'test' in composition and composition['test']:
                sorted_residues = sorted(
                    [[aa, count] for aa, count in composition['test'].items()],
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
            'unique_aas',
            'has_diff_majority',
            'aa_diversity',
            'majority_residue',
            'composition'
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
            
            # Use sequence_utils to process the sequence data for this site
            # Convert composition dict of Counters to the format expected by sequence_utils
            group_sequences = {}
            for group_name, aa_counter in composition.items():
                if aa_counter:  # Only process non-empty counters
                    group_sequences[group_name] = aa_counter
            
            # Process sequence data using our improved function
            if group_sequences:
                site_metrics = su.process_sequence_data(group_sequences)
                
                # Add site metrics to comparison data for each group
                for group_name in self._comparison_groups:
                    if group_name in group_sequences:
                        # Initialize group data if not exists
                        if group_name not in site_comparison_data:
                            site_comparison_data[group_name] = {}
                            
                        # Add metrics for this group
                        site_comparison_data[group_name].update({
                            'unique_aas': ','.join(site_metrics.get('unique_aas', {}).get(group_name, [])) or 'NA',
                            'has_diff_majority': site_metrics.get('has_diff_majority', False),
                            'aa_diversity': site_metrics.get(f'{group_name}_diversity', 0),
                            'majority_residue': site_metrics.get(f'{group_name}_majority', '-')
                        })
                        
                        # Add composition data
                        if f'{group_name}_composition' in site_metrics:
                            formatted_comp = su.format_composition(site_metrics[f'{group_name}_composition'])
                            site_comparison_data[group_name]['composition'] = formatted_comp
            
            # Add comparison data for any groups that weren't processed
            for group in self._comparison_groups:
                if group not in site_comparison_data and group in composition and composition[group]:
                    site_comparison_data[group] = {
                        'unique_aas': 'NA',
                        'has_diff_majority': False,
                        'aa_diversity': 0,
                        'majority_residue': '-',
                        'composition': ''
                    }
                    
            # Only add sites with comparison data
            if site_comparison_data:
                comparison_data[site_num] = site_comparison_data
        
        return comparison_data
