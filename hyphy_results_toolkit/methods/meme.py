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
        if 'MLE' in results and 'content' in results['MLE']:
            for row in results['MLE']['content']['0']:
                p_value = float(row[7])  # P-value for episodic selection
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
        
        if 'MLE' in results and 'content' in results['MLE']:
            for row in results['MLE']['content']['0']:
                site = int(row[0])
                site_results[site] = {
                    'meme_alpha': float(row[2]),     # Synonymous rate
                    'meme_beta_neg': float(row[3]),  # Non-synonymous rate for negative selection
                    'meme_beta_plus': float(row[4]), # Non-synonymous rate for positive selection
                    'meme_weight': float(row[5]),    # Weight of positive selection class
                    'meme_pvalue': float(row[7]),    # P-value for episodic selection
                    'meme_selection': 'episodic' if float(row[7]) <= 0.05 else 'none'
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
    def get_site_fields(clades: List[str] = None) -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'meme_alpha',
            'meme_beta_neg',
            'meme_beta_plus',
            'meme_weight',
            'meme_pvalue',
            'meme_selection'
        ]
