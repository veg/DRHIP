"""
PRIME (PRoperty Informed Models of Evolution) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class PrimeMethod(HyPhyMethod):
    """Implementation of PRIME analysis processing."""
    
    PROPERTIES = [
        'Chemical composition',
        'Coil propensity',
        'Compressibility',
        'Hydropathy',
        'Molecular volume',
        'Polarity',
        'Refractivity',
        'Surface composition'
    ]
    
    def __init__(self):
        """Initialize PRIME method."""
        super().__init__("PRIME", "PRIME.json")
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process PRIME results.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        processed = {}
        
        # Process sites under selection for each property
        for prop in self.PROPERTIES:
            sites_conserved = 0
            sites_altered = 0
            
            if self.has_mle_content(results) and self.has_mle_headers(results):
                # Get header indices
                header_indices = self.get_header_indices(results)
                
                # Get index for p-value
                pvalue_index = self.get_column_index(header_indices, 'p-value', 9)
                
                for row in results['MLE']['content']['0']:
                    p_value = float(row[pvalue_index])
                    if p_value <= 0.05:
                        sites_conserved += 1
                    else:
                        sites_altered += 1
            
            # We don't need to store property-specific counts in the summary
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific PRIME data.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)
            
            # Get indices for the values we need
            site_index = self.get_column_index(header_indices, 'Site', 0)
            pvalue_index = self.get_column_index(header_indices, 'p-value', 9)
            
            for row in results['MLE']['content']['0']:
                site = int(row[site_index]) if site_index >= 0 else 0
                if site_index < 0:
                    site += 1  # If site index not found, increment counter
                
                p_value = float(row[pvalue_index])
                
                site_results[site] = {
                    # Only include the marker field
                    'prime_marker': 'overall' if p_value <= 0.05 else '-'
                }
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []  # No summary fields needed
    
    @staticmethod
    def get_site_fields(comparison_groups: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'prime_marker'
        ]
