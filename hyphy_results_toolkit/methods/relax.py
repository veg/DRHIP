"""
RELAX (Relaxation) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class RelaxMethod(HyPhyMethod):
    """Implementation of RELAX analysis processing."""
    
    def __init__(self):
        """Initialize RELAX method."""
        super().__init__("RELAX", "RELAX.json")
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process RELAX results.
        
        Args:
            results: Raw RELAX results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        test_results = results['test results']
        return {
            'relax_pvalue': test_results['p-value'],
            'relax_k': test_results['relaxation or intensification parameter'],
            'relax_lrt': test_results['LRT'],
            'relax_k_significant': test_results['p-value'] <= 0.05
        }
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'relax_pvalue',
            'relax_k',
            'relax_lrt',
            'relax_k_significant'
        ]
