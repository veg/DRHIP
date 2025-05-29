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
        processed = {
            'prime_sites_tested': len(results.get('MLE', {}).get('content', {}).get('0', [])),
            'prime_version': results.get('version', ''),
            'prime_timestamp': results.get('timestamp', '')
        }
        
        # Process sites under selection for each property
        for prop in self.PROPERTIES:
            sites_conserved = 0
            sites_altered = 0
            
            if 'MLE' in results and 'content' in results['MLE']:
                for row in results['MLE']['content']['0']:
                    p_value = row[9]
                    if p_value <= 0.05:
                        sites_conserved += 1
                    else:
                        sites_altered += 1
            
            prop_key = prop.lower().replace(' ', '_')
            processed[f'prime_{prop_key}_conserved'] = sites_conserved
            processed[f'prime_{prop_key}_altered'] = sites_altered
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific PRIME data.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        
        if 'MLE' in results and 'content' in results['MLE']:
            site_index = 1
            for row in results['MLE']['content']['0']:
                site_dict = {}
                
                for prop in self.PROPERTIES:
                    p_value = row[9]
                    
                    site_dict.update({
                        f'prime_{site_index}_pvalue': p_value,
                        f'prime_{site_index}_conserved': p_value <= 0.05,
                        f'prime_{site_index}_altered': p_value > 0.05
                    })
                
                site_results[site_index-1] = site_dict
                site_index += 1
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        fields = [
            'prime_sites_tested',
            'prime_version',
            'prime_timestamp'
        ]
        
        # Add fields for each property
        for prop in PrimeMethod.PROPERTIES:
            prop_key = prop.lower().replace(' ', '_')
            fields.extend([
                f'prime_{prop_key}_conserved',
                f'prime_{prop_key}_altered'
            ])
        
        return fields
    
    @staticmethod
    def get_site_fields(comparison_groups: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        # For PRIME, we'll return an empty list since we handle site fields differently
        # The actual site fields will be determined at runtime based on the data
        return []
