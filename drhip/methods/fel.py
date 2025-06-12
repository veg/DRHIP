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
        
        # Calculate selection statistics using the helper
        selection_stats = self.calculate_selection_counts(
            results,
            alpha_col='alpha',
            beta_col='beta',
            pvalue_col='p-value',
            significance=0.05
        )
        
        # Update processed results with selection statistics
        processed.update(selection_stats)
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
        """Process site-specific FEL data.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        # Define column names mapping
        column_names = {
            'alpha': 'alpha',
            'beta': 'beta',
            'p-value': 'p-value'
        }
        
        # Define a function to process each row
        def process_row(site_idx, row, column_indices):
            alpha_index = column_indices['alpha']
            beta_index = column_indices['beta']
            pvalue_index = column_indices['p-value']
            
            # Check if we have valid indices and data
            if (alpha_index < 0 or beta_index < 0 or pvalue_index < 0 or
                alpha_index >= len(row) or beta_index >= len(row) or pvalue_index >= len(row)):
                # Return NA values for missing data
                return {
                    'fel_selection': 'NA'
                }
            
            try:
                alpha = float(row[alpha_index])
                beta = float(row[beta_index])
                pvalue = float(row[pvalue_index])
                
                # Determine selection type
                selection_type = (
                    'positive' if beta > alpha and pvalue <= 0.05
                    else 'negative' if beta < alpha and pvalue <= 0.05
                    else 'neutral'
                )
                
                # Return site data - only include site-specific fields
                return {
                    'fel_selection': selection_type if pvalue <= 0.05 else 'neutral'
                }
            except (ValueError, TypeError):
                # Return NA values for malformed data
                return {
                    'fel_selection': 'NA'
                }
        
        # Use the helper to process site data
        return self.process_site_mle_data(results, column_names, process_row)
    
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
        
    @staticmethod
    def get_comparison_group_fields() -> List[str]:
        """Get list of fields that are specific to comparison groups."""
        return []  # FEL doesn't have comparison group-specific fields
        
    def process_comparison_site_data(self, results: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from FEL results.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # FEL doesn't have comparison group-specific site data
        return {}
