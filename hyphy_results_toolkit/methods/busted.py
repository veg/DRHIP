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
            'busted_pvalue': test_results['p-value'],
            'busted_lrt': test_results['LRT'],
            'busted_evidence': test_results['p-value'] <= 0.05
        }
        
        # Get omega3 distribution
        omega3 = self.get_omega3(results)
        processed['busted_omega3'] = omega3['omega']
        # Use 'proportion' key for CAPHEINE format or 'weight' for original format
        if 'weight' in omega3:
            processed['busted_omega3_weight'] = omega3['weight']
        elif 'proportion' in omega3:
            processed['busted_omega3_weight'] = omega3['proportion']
        else:
            processed['busted_omega3_weight'] = 0.0
        
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
            
            processed['total_branch_length'] = total_branch_length
            
            # Calculate tree length ratio if reference exists
            if 'reference' in results.get('branch attributes', {}):
                ref_branch_lengths = results['branch attributes']['reference']
                ref_total_branch_length = 0
                
                # Check if the reference branches have a 'length' key
                if all('length' in branch for branch in ref_branch_lengths.values()):
                    ref_total_branch_length = sum(float(branch['length']) for branch in ref_branch_lengths.values())
                # Otherwise, try to use MG94xREV values as branch lengths
                else:
                    for node_name, node_data in ref_branch_lengths.items():
                        if 'MG94xREV with separate rates for branch sets' in node_data:
                            ref_total_branch_length += float(node_data['MG94xREV with separate rates for branch sets'])
                
                if ref_total_branch_length > 0:
                    processed['tree_length_ratio'] = total_branch_length / ref_total_branch_length
                else:
                    processed['tree_length_ratio'] = 1.0
            else:
                processed['tree_length_ratio'] = 1.0
        except Exception:
            # If anything goes wrong with branch length calculations, use default values
            processed['total_branch_length'] = 0.0
            processed['tree_length_ratio'] = 1.0
        
        return processed
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'busted_pvalue',
            'busted_lrt',
            'busted_evidence',
            'busted_omega3',
            'busted_omega3_weight',
            'total_branch_length',
            'tree_length_ratio'
        ]
