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
    
    def validate_input_json(self, results: Dict[str, Any]) -> List[str]:
        missing = []
        missing.extend(self.validate_required_paths(results, ['MLE.headers', 'MLE.content.0']))
        # If headers exist, ensure key columns are present
        missing_headers = self.missing_mle_headers(results, ['alpha', 'beta', 'p-value'])
        if missing_headers:
            missing.extend([f"MLE.headers:{h}" for h in missing_headers])
        return missing
        
    def calculate_negative_sites(self, results: Dict[str, Any], 
                                alpha_col: str = 'alpha', 
                                beta_col: str = 'beta', 
                                pvalue_col: str = 'p-value',
                                significance: float = 0.05) -> Dict[str, int]:
        """Calculate negative selection statistics from MLE data.
        
        This helper calculates FEL negative selection statistics:
        - negative_sites: Number of sites under negative selection (beta < alpha, p <= significance)
        
        Args:
            results: Raw results dictionary from JSON file
            alpha_col: Name of the column containing alpha (synonymous rate) values
            beta_col: Name of the column containing beta (non-synonymous rate) values
            pvalue_col: Name of the column containing p-values
            significance: P-value threshold for significance (default: 0.05)
            
        Returns:
            Dictionary with negative selection statistics
        """
        stats = {
            'negative_sites': 0
        }
        
        # Check if required data is available
        if not self.has_mle_content(results) or not self.has_mle_headers(results):
            return stats
            
        # Get header indices
        header_indices = self.get_header_indices(results)
        
        # Get indices for the values we need
        alpha_index = self.get_column_index(header_indices, alpha_col, -1)
        beta_index = self.get_column_index(header_indices, beta_col, -1)
        pvalue_index = self.get_column_index(header_indices, pvalue_col, -1)
        
        # Check if we have valid column indices
        if alpha_index < 0 or beta_index < 0 or pvalue_index < 0:
            return stats
        
        # Process each site
        for row in results['MLE']['content']['0']:
            try:
                # Extract values
                alpha = float(row[alpha_index])    # Alpha (synonymous rate)
                beta = float(row[beta_index])      # Beta (non-synonymous rate)
                p_value = float(row[pvalue_index]) # P-value
                
                # Count sites based on selection criteria
                if p_value <= significance and beta < alpha:
                    stats['negative_sites'] += 1
            except (ValueError, IndexError, TypeError):
                # Skip this site if there's an error
                continue
        
        return stats
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process FEL results.
        
        Args:
            results: Raw FEL results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        # Extract common fields (N, T, sites)
        processed = self.extract_common_fields(results)
        
        # Calculate negative selection statistics using the helper
        selection_stats = self.calculate_negative_sites(
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
