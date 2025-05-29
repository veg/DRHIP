"""
FEL (Fixed Effects Likelihood) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class FelMethod(HyPhyMethod):
    """Implementation of FEL analysis processing."""
    
    def __init__(self):
        """Initialize FEL method."""
        super().__init__("FEL", "FEL.json")
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process FEL results.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        processed = {
            'fel_sites_tested': len(results.get('MLE', {}).get('content', {}).get('0', [])),
            'fel_version': results.get('version', ''),
            'fel_timestamp': results.get('timestamp', '')
        }
        
        # Process tested sites
        sites_under_selection = 0
        sites_under_negative_selection = 0
        
        if 'MLE' in results and 'content' in results['MLE']:
            for row in results['MLE']['content']['0']:
                alpha = float(row[0])  # Alpha (synonymous rate)
                beta = float(row[1])   # Beta (non-synonymous rate)
                p_value = float(row[4])  # P-value
                
                if p_value <= 0.05:
                    if beta > alpha:
                        sites_under_selection += 1
                    elif beta < alpha:
                        sites_under_negative_selection += 1
        
        processed.update({
            'fel_sites_positive_selection': sites_under_selection,
            'fel_sites_negative_selection': sites_under_negative_selection
        })
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific FEL data.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        
        if 'MLE' in results and 'content' in results['MLE']:
            for row in results['MLE']['content']['0']:
                site = int(row[0])
                site_results[site] = {
                    'fel_alpha': float(row[0]),      # Synonymous rate
                    'fel_beta': float(row[1]),       # Non-synonymous rate
                    'fel_pvalue': float(row[4]),     # P-value
                    'fel_selection': (
                        'positive' if float(row[1]) > float(row[0]) and float(row[4]) <= 0.05
                        else 'negative' if float(row[1]) < float(row[0]) and float(row[4]) <= 0.05
                        else 'neutral'
                    )
                }
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'fel_sites_tested',
            'fel_sites_positive_selection',
            'fel_sites_negative_selection',
            'fel_version',
            'fel_timestamp'
        ]
    
    @staticmethod
    def get_site_fields(comparison_groups: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'fel_alpha',
            'fel_beta',
            'fel_pvalue',
            'fel_selection'
        ]
