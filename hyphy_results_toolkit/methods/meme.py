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
        processed = {
            'meme_sites_tested': len(results.get('MLE', {}).get('content', {}).get('0', [])),
            'meme_version': results.get('version', ''),
            'meme_timestamp': results.get('timestamp', '')
        }
        
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
        
        processed['meme_sites_selection'] = sites_under_selection
        
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
            alpha_index = self.get_column_index(header_indices, 'alpha', 2)
            beta_neg_index = self.get_column_index(header_indices, 'beta-', 3)
            beta_plus_index = self.get_column_index(header_indices, 'beta+', 4)
            weight_index = self.get_column_index(header_indices, 'weight', 5)
            pvalue_index = self.get_column_index(header_indices, 'p-value', 7)
            
            for row in results['MLE']['content']['0']:
                site = int(row[site_index])
                alpha = float(row[alpha_index])
                beta_neg = float(row[beta_neg_index])
                beta_plus = float(row[beta_plus_index])
                weight = float(row[weight_index])
                pvalue = float(row[pvalue_index])
                
                site_results[site] = {
                    'meme_alpha': alpha,     # Synonymous rate
                    'meme_beta_neg': beta_neg,  # Non-synonymous rate for negative selection
                    'meme_beta_plus': beta_plus, # Non-synonymous rate for positive selection
                    'meme_weight': weight,    # Weight of positive selection class
                    'meme_pvalue': pvalue,    # P-value for episodic selection
                    'meme_selection': 'episodic' if pvalue <= 0.05 else 'none'
                }
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'meme_sites_tested',
            'meme_sites_selection',
            'meme_version',
            'meme_timestamp'
        ]
    
    @staticmethod
    def get_site_fields(comparison_groups: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'meme_alpha',
            'meme_beta_neg',
            'meme_beta_plus',
            'meme_weight',
            'meme_pvalue',
            'meme_selection'
        ]
