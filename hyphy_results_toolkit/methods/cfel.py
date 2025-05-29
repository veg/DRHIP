"""
CFEL (Contrast-FEL) method implementation.
"""

from typing import Dict, Any, List, Tuple

from .base import HyPhyMethod

class CfelMethod(HyPhyMethod):
    """Implementation of CFEL analysis processing."""
    
    def __init__(self):
        """Initialize CFEL method."""
        super().__init__("CFEL", "CFEL.json")
        self._header_map = {}
        self._beta_idx_map = {}
        self._subs_idx_map = {}
    
    def _build_column_maps(self, headers: List[Tuple[str, str]], comparison_groups: List[str]) -> None:
        """Build column lookup maps for the CFEL data.
        
        Args:
            headers: List of (short_label, description) tuples
            comparison_groups: List of comparison group names
        """
        self._header_map = {short_label: i for i, (short_label, _) in enumerate(headers)}
        self._beta_idx_map = {g: self._header_map[f"beta ({g})"] for g in comparison_groups}
        self._subs_idx_map = {g: self._header_map[f"subs ({g})"] for g in comparison_groups}
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process CFEL results.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        # Get the tag for each branch
        tested = results['tested']['0']
        by_type: Dict[str, List[str]] = {}
        
        # Group branch names by tag
        for branch, tag in tested.items():
            if tag not in by_type:
                by_type[tag] = []
            by_type[tag].append(branch)
        
        # Build column lookup maps
        headers = results["MLE"]["headers"]
        self._build_column_maps(headers, list(by_type.keys()))
        
        return {
            'by_type': by_type,
            'data_rows': results["MLE"]["content"]["0"],
            'beta_idx_map': self._beta_idx_map,
            'subs_idx_map': self._subs_idx_map
        }
    
    def process_site_data(self, site_data: Dict[str, Any], comparison_groups: List[str]) -> Dict[str, Any]:
        """Process site-specific CFEL data.
        
        Args:
            site_data: Site data from process_results
            comparison_groups: List of comparison groups to process
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        data_rows = site_data['data_rows']
        beta_idx_map = site_data['beta_idx_map']
        subs_idx_map = site_data['subs_idx_map']
        
        for row in data_rows:
            site = int(row[0])  # First column is the site index
            site_dict = {}
            
            # Extract beta values and substitution counts for each comparison group
            for group in comparison_groups:
                site_dict[f'beta_{group}'] = float(row[beta_idx_map[group]])
                site_dict[f'subs_{group}'] = float(row[subs_idx_map[group]])
            
            site_results[site] = site_dict
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'nt_conserved',
            'aa_conserved'
        ]
    
    @staticmethod
    def get_site_fields(comparison_groups: List[str]) -> List[str]:
        """Get list of site-specific fields produced by this method.
        
        Args:
            comparison_groups: List of comparison groups to generate fields for
            
        Returns:
            List of field names
        """
        fields = []
        for group in comparison_groups:
            fields.extend([
                f'beta_{group}',
                f'subs_{group}'
            ])
        return fields
