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
            # Get the headers directly from the results
            headers = results['MLE']['headers']
            
            # Extract property tags from headers, similar to End2End-DENV
            prime_tags = []
            p_indices = []
            
            for i, (short_label, description) in enumerate(headers):
                if short_label == 'p-value':
                    p_indices.append(i)
                    prime_tags.append((0, 'overall'))
                elif short_label.startswith('p') and short_label[1:].isdigit():
                    p_indices.append(i)
                    # Extract property name from description
                    property_name = description.split(' ')[-1] if ' ' in description else short_label
                    prime_tags.append((int(short_label[1:]), property_name))
            
            # Sort tags by their index
            prime_tags = [tag for _, tag in sorted(prime_tags, key=lambda x: x[0])]
            
            # Process each site
            for row_idx, row in enumerate(results['MLE']['content']['0']):
                site = row_idx  # Use row index as site number
                
                # Get p-values for all properties
                pvals = [float(row[idx]) for idx in p_indices]
                
                # If any p-value is significant
                if min(pvals) <= 0.05:
                    # Collect significant property tags
                    significant_tags = [prime_tags[i] for i, pval in enumerate(pvals) if pval <= 0.05]
                    
                    site_results[site] = {
                        'prime_marker': ','.join(significant_tags) if significant_tags else '-'
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
            'prime_marker'
        ]
