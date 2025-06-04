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
        # Extract common fields (N, T, sites)
        processed = self.extract_common_fields(results)
        
        # Process tested sites
        sites_under_selection = 0
        sites_under_negative_selection = 0
        diff_sites = 0
        
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)
            
            # Get indices for the values we need
            alpha_index = self.get_column_index(header_indices, 'alpha', 0)
            beta_index = self.get_column_index(header_indices, 'beta', 1)
            pvalue_index = self.get_column_index(header_indices, 'p-value', 4)
            
            for row in results['MLE']['content']['0']:
                alpha = float(row[alpha_index])  # Alpha (synonymous rate)
                beta = float(row[beta_index])    # Beta (non-synonymous rate)
                p_value = float(row[pvalue_index])  # P-value
                
                if p_value <= 0.05:
                    diff_sites += 1
                    if beta > alpha:
                        sites_under_selection += 1
                    elif beta < alpha:
                        sites_under_negative_selection += 1
        
        processed.update({
            'positive_sites': sites_under_selection,  
            'negative_sites': sites_under_negative_selection,
            'diff_sites': diff_sites
        })
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
        """Process site-specific FEL data.
        
        Args:
            results: Raw FEL results dictionary
            
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
        
        # Get indices for the values we need
        alpha_index = self.get_column_index(header_indices, 'alpha', 0)
        beta_index = self.get_column_index(header_indices, 'beta', 1)
        pvalue_index = self.get_column_index(header_indices, 'p-value', 4)
        
        # Process each site (row index + 1 is the site number)
        for site_idx, row in enumerate(results['MLE']['content']['0'], 1):
            # Check if we have valid data for this site
            has_valid_data = True
            
            # Ensure we have valid values
            try:
                if alpha_index < len(row) and beta_index < len(row) and pvalue_index < len(row):
                    alpha = float(row[alpha_index])
                    beta = float(row[beta_index])
                    pvalue = float(row[pvalue_index])
                else:
                    has_valid_data = False
            except (ValueError, TypeError):
                has_valid_data = False
            
            if not has_valid_data:
                # Use NA for missing or malformed data
                site_data = {
                    'fel_selection': 'NA',
                    'intensified_positive_selection': 'NA',
                    'cfel_marker': 'NA'
                }
            else:
                # Determine selection type
                selection_type = (
                    'positive' if beta > alpha and pvalue <= 0.05
                    else 'negative' if beta < alpha and pvalue <= 0.05
                    else 'neutral'
                )
                
                # Initialize site data
                site_data = {
                    'fel_selection': selection_type if pvalue <= 0.05 else 'neutral',
                    'intensified_positive_selection': beta > alpha and pvalue <= 0.05,
                    'cfel_marker': f"{pvalue:.3f}" if pvalue <= 0.05 else "-"
                }
            
            site_results[site_idx] = site_data
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'positive_sites',
            'negative_sites'
        ]
    
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'fel_selection',
            'intensified_positive_selection',
            'cfel_marker'
        ]
