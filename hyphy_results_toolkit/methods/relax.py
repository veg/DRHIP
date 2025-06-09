"""
RELAX (Relaxation) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class RelaxMethod(HyPhyMethod):
    """Implementation of RELAX analysis processing."""
    
    def __init__(self):
        """Initialize RELAX method."""
        super().__init__("RELAX", "RELAX.json")
        self._comparison_groups = None
    
    def process_comparison_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process RELAX results.
        
        Args:
            results: Raw RELAX results dictionary
            
        Returns:
            Processed results with standardized keys
            
        Raises:
            ValueError: If comparison groups are not set
        """
        # Ensure comparison groups are set
        if not self._comparison_groups:
            raise ValueError("Comparison groups must be set before processing RELAX results")
            
        # Check if results contain the required data
        if not results or 'test results' not in results:
            # Return NA for all fields if data is missing
            processed = {
                'RELAX_overall_pval': 'NA',
                'RELAX_K': 'NA'
            }
            
            # Add NA for all group-specific fields
            for group in self._comparison_groups:
                processed[f'RELAX_K_{group}'] = 'NA'
                
            return processed
        
        test_results = results['test results']
        processed = {
            'RELAX_overall_pval': test_results['p-value'],  # Overall p-value field
            'RELAX_K': test_results['relaxation or intensification parameter']  # Overall K parameter
        }
        
        # Extract group-specific K parameter values if available
        if isinstance(test_results['relaxation or intensification parameter'], dict):
            # For each comparison group, extract its K parameter value
            for group in self._comparison_groups:
                if group in test_results['relaxation or intensification parameter']:
                    processed[f'RELAX_K_{group}'] = test_results['relaxation or intensification parameter'][group]
        
        return processed
    
    def get_comparison_group_summary_fields(self) -> List[str]:
        """Get list of summary fields produced by this method.
        
        Raises:
            ValueError: If comparison groups are not set
        """
        # Ensure comparison groups are set
        if not self._comparison_groups:
            raise ValueError("Comparison groups must be set before getting summary fields")
            
        base_fields = [
            'RELAX_overall_pval',
            'RELAX_K'
        ]
        
        # Add group-specific fields
        for group in self._comparison_groups:
            base_fields.append(f'RELAX_K_{group}')
            
        return base_fields
        
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return []  # RELAX doesn't have site-specific fields

    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []  # RELAX doesn't have summary fields
        
    @staticmethod
    def get_comparison_group_site_fields() -> List[str]:
        """Get list of fields that are specific to comparison groups."""
        return []  # RELAX doesn't have comparison group-specific site fields
        
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process RELAX results.
        
        Args:
            results: Raw RELAX results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        # RELAX doesn't use standard summary fields, only comparison group fields
        return {}

    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific RELAX data.
        
        Args:
            results: Raw RELAX results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        # RELAX doesn't have site-specific data
        return {}
        
    def process_comparison_site_data(self, results: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from RELAX results.
        
        Args:
            results: Raw RELAX results dictionary
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # RELAX doesn't have comparison group-specific site data
        return {}
