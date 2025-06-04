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
        
        # Check if required data is available
        if not self.has_mle_content(results) or not self.has_mle_headers(results):
            # If critical data is missing, return empty results
            return site_results
            
        # Get header indices
        header_indices = self.get_header_indices(results)
        
        # Get indices for the values we need (with defaults matching original hardcoded indices)
        site_index = self.get_column_index(header_indices, 'Site', 0)
        pvalue_index = self.get_column_index(header_indices, 'p-value', 7)
        
        for row in results['MLE']['content']['0']:
            # Get site index with error handling
            try:
                site = int(row[site_index]) if row and len(row) > site_index else 0
            except (ValueError, TypeError, IndexError):
                continue  # Skip this row if we can't get a valid site index
            
            # Check if we have valid data for this site
            has_valid_data = True
            
            # Ensure we have a valid p-value
            try:
                if pvalue_index < len(row):
                    pvalue = float(row[pvalue_index])
                else:
                    has_valid_data = False
            except (ValueError, TypeError):
                has_valid_data = False
            
            # Set the marker based on data validity and significance
            if not has_valid_data:
                marker = "NA"  # Use NA for missing or malformed data
            elif pvalue <= 0.05:
                marker = f'{pvalue:.3f}'  # Format p-value for significant sites
            else:
                marker = '-'  # Use dash for non-significant sites
            
            site_results[site] = {
                'meme_marker': marker
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
