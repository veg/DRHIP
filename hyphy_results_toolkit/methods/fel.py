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
        processed = {}
        
        # Extract number of sites from input data
        if 'input' in results and 'number of sites' in results['input']:
            processed['sites'] = results['input']['number of sites']
        elif 'input' in results and 'sites' in results['input']:
            processed['sites'] = results['input']['sites']
        elif self.has_mle_content(results):
            # If not explicitly provided, use the number of rows in the MLE content
            processed['sites'] = len(results['MLE']['content']['0'])
        else:
            processed['sites'] = 0
        
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
        
        # Extract sequence count if available
        if 'input' in results and 'number of sequences' in results['input']:
            processed['N'] = results['input']['number of sequences']
        
        # Extract branch length information if available
        if 'branch attributes' in results and '0' in results['branch attributes']:
            try:
                branch_lengths = [float(branch.get('length', 0)) for branch in results['branch attributes']['0'].values()]
                processed['T'] = sum(branch_lengths)
            except (ValueError, TypeError, KeyError):
                pass
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
        """Process site-specific FEL data.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)
            
            # Get indices for the values we need
            alpha_index = self.get_column_index(header_indices, 'alpha', 0)
            beta_index = self.get_column_index(header_indices, 'beta', 1)
            pvalue_index = self.get_column_index(header_indices, 'p-value', 4)
            
            # Process each site (row index + 1 is the site number)
            for site_idx, row in enumerate(results['MLE']['content']['0'], 1):
                alpha = float(row[alpha_index])
                beta = float(row[beta_index])
                pvalue = float(row[pvalue_index])
                
                # Determine selection type
                selection_type = (
                    'positive' if beta > alpha and pvalue <= 0.05
                    else 'negative' if beta < alpha and pvalue <= 0.05
                    else 'neutral'
                )
                
                # In End2End-DENV, FEL stores [p-value, is_positive]
                # where is_positive is True if beta > alpha (positive selection)
                if pvalue <= 0.05:
                    site_results[site_idx] = {
                        'fel_selection': selection_type
                    }
                    
                    # Update site-specific fields for the desired output format
                    if selection_type == 'positive':
                        # Mark this site as under positive selection
                        site_results[site_idx]['intensified_positive_selection'] = True
        
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
            'fel_selection'
        ]
