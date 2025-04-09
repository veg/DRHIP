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
                    prop_data = row[1].get(prop, {})
                    p_value = float(prop_data.get('p-value', 1.0))
                    if p_value <= 0.05:
                        if prop_data.get('preserved', False):
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
            for row in results['MLE']['content']['0']:
                site = int(row[0])
                site_dict = {}
                
                for prop in self.PROPERTIES:
                    prop_data = row[1].get(prop, {})
                    prop_key = prop.lower().replace(' ', '_')
                    
                    site_dict.update({
                        f'prime_{prop_key}_pvalue': float(prop_data.get('p-value', 1.0)),
                        f'prime_{prop_key}_preserved': prop_data.get('preserved', False),
                        f'prime_{prop_key}_altered': not prop_data.get('preserved', True)
                    })
                
                site_results[site] = site_dict
        
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
    def get_site_fields(clades: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        fields = []
        
        # Add fields for each property
        for prop in PrimeMethod.PROPERTIES:
            prop_key = prop.lower().replace(' ', '_')
            fields.extend([
                f'prime_{prop_key}_pvalue',
                f'prime_{prop_key}_preserved',
                f'prime_{prop_key}_altered'
            ])
        
        return fields
