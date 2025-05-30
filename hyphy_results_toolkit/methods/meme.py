"""
MEME (Mixed Effects Model of Evolution) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class MemeMethod(HyPhyMethod):
    """Implementation of MEME analysis processing."""
    
    def __init__(self):
        """Initialize MEME method."""
        super().__init__("MEME", "MEME.json")
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process MEME results.
        
        Args:
            results: Raw MEME results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        processed = {}
        
        # Process sites under selection
        sites_under_selection = 0
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)
            
            # Get index for p-value
            pvalue_index = self.get_column_index(header_indices, 'p-value', 7)  # P-value for episodic selection
            
            for row in results['MLE']['content']['0']:
                p_value = float(row[pvalue_index])
                if p_value <= 0.05:
                    sites_under_selection += 1
        
        # We don't need to store this in the summary
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific MEME data.
        
        Args:
            results: Raw MEME results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)
            
            # Get indices for the values we need (with defaults matching original hardcoded indices)
            site_index = self.get_column_index(header_indices, 'Site', 0)
            pvalue_index = self.get_column_index(header_indices, 'p-value', 7)
            
            for row in results['MLE']['content']['0']:
                site = int(row[site_index])
                pvalue = float(row[pvalue_index])
                
                site_results[site] = {
                    'meme_marker': f'{pvalue:.3f}' if pvalue <= 0.05 else '-'
                }
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []  # No summary fields needed
    
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'meme_marker'
        ]
