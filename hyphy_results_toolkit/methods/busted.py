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
        processed['busted_omega3_weight'] = omega3['weight']
        
        # Get branch lengths
        branch_lengths = results['branch attributes']['0']
        total_branch_length = sum(float(branch['length']) for branch in branch_lengths.values())
        processed['total_branch_length'] = total_branch_length
        
        # Calculate tree length ratio
        if 'reference' in results['branch attributes']:
            ref_branch_lengths = results['branch attributes']['reference']
            ref_total_length = sum(float(branch['length']) for branch in ref_branch_lengths.values())
            processed['tree_length_ratio'] = total_branch_length / ref_total_length
        else:
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
