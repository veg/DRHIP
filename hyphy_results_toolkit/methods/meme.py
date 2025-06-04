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
        # Define column names mapping
        column_names = {
            'site': 'Site',
            'p-value': 'p-value'
        }
        
        # Define a function to process each row
        def process_row(site_idx, row, column_indices):
            pvalue_index = column_indices['p-value']
            
            # Check if we have valid indices and data
            if pvalue_index < 0 or pvalue_index >= len(row):
                return {
                    'meme_marker': 'NA'  # Use NA for missing data
                }
            
            try:
                pvalue = float(row[pvalue_index])
                
                # Set the marker based on significance
                if pvalue <= 0.05:
                    marker = f'{pvalue:.3f}'  # Format p-value for significant sites
                else:
                    marker = '-'  # Use dash for non-significant sites
                
                return {
                    'meme_marker': marker
                }
            except (ValueError, TypeError):
                # Return NA for malformed data
                return {
                    'meme_marker': 'NA'
                }
        
        # Use the helper to process site data
        return self.process_site_mle_data(results, column_names, process_row)
    
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
