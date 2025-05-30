"""
BUSTED (Branch-Site Unrestricted Statistical Test for Episodic Diversification) method implementation.
"""

from typing import Dict, Any, List

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
        test_results = results['test results']
        processed = {
            'BUSTED_pval': test_results['p-value']  # Only keep the p-value field
        }
        
        # Get omega3 distribution
        omega3 = self.get_omega3(results)
        processed['BUSTED_omega3'] = omega3['omega']  # Renamed to match desired format
        
        # Get proportion of sites in omega3 category
        if 'weight' in omega3:
            processed['BUSTED_prop_sites_in_omega3'] = omega3['weight']  # Renamed to match desired format
        elif 'proportion' in omega3:
            processed['BUSTED_prop_sites_in_omega3'] = omega3['proportion']  # Renamed to match desired format
        else:
            processed['BUSTED_prop_sites_in_omega3'] = 0.0
        
        # Get branch lengths - handle different formats
        try:
            # Try to calculate branch lengths
            total_branch_length = 0
            
            # Check if we have the original format with 'length' key
            if '0' in results.get('branch attributes', {}):
                branch_lengths = results['branch attributes']['0']
                # Check if the branches have a 'length' key
                if all('length' in branch for branch in branch_lengths.values()):
                    total_branch_length = sum(float(branch['length']) for branch in branch_lengths.values())
                # Otherwise, try to use MG94xREV values as branch lengths
                else:
                    for node_name, node_data in branch_lengths.items():
                        if 'MG94xREV with separate rates for branch sets' in node_data:
                            total_branch_length += float(node_data['MG94xREV with separate rates for branch sets'])
            
            processed['T'] = total_branch_length  # Use T for total branch length
            
            # Count number of sequences from tree
            try:
                # Try to get number of sequences from the tree
                if 'tested' in results and 'sequences' in results['tested']:
                    processed['N'] = len(results['tested']['sequences'])
                else:
                    # Estimate number of sequences from branch attributes
                    # This is a rough estimate - half the number of branches in a bifurcating tree
                    processed['N'] = len(branch_lengths) // 2
            except Exception:
                processed['N'] = 0
            
            # Get number of sites
            if 'input' in results and 'sites' in results['input']:
                processed['sites'] = results['input']['sites']
            else:
                processed['sites'] = 0
                
            # Calculate dN/dS if possible
            if 'fits' in results and 'Unconstrained model' in results['fits']:
                model_fit = results['fits']['Unconstrained model']
                if 'Rate Distributions' in model_fit and 'global' in model_fit['Rate Distributions']:
                    rates = model_fit['Rate Distributions']['global']
                    # Calculate weighted average of omega values
                    omega_sum = 0.0
                    for rate in rates:
                        omega_sum += rate['omega'] * rate.get('weight', rate.get('proportion', 0.0))
                    processed['dN/dS'] = omega_sum
                else:
                    processed['dN/dS'] = 0.0
            else:
                processed['dN/dS'] = 0.0
            
        except Exception as e:
            # If anything goes wrong with calculations, use default values
            print(f"Error processing BUSTED results: {e}")
            processed['T'] = 0.0
            processed['N'] = 0
            processed['sites'] = 0
            processed['dN/dS'] = 0.0
        
        return processed
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'BUSTED_pval',
            'BUSTED_omega3',
            'BUSTED_prop_sites_in_omega3'
        ]
