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
        processed = {
            'RELAX_overall_pval': test_results['p-value'],  # p-value field
            'RELAX_K': test_results['relaxation or intensification parameter']  # K parameter
        }
        
        return processed
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'RELAX_overall_pval',
            'RELAX_K'
        ]
