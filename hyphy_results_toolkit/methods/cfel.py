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
    
    def _build_column_maps(self, headers: List[Tuple[str, str]], clades: List[str]) -> None:
        """Build column lookup maps for the CFEL data.
        
        Args:
            headers: List of (short_label, description) tuples
            clades: List of clade names
        """
        self._header_map = {short_label: i for i, (short_label, _) in enumerate(headers)}
        self._beta_idx_map = {c: self._header_map[f"beta ({c})"] for c in clades}
        self._subs_idx_map = {c: self._header_map[f"subs ({c})"] for c in clades}
    
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
    
    def process_site_data(self, site_data: Dict[str, Any], clades: List[str]) -> Dict[str, Any]:
        """Process site-specific CFEL data.
        
        Args:
            site_data: Site data from process_results
            clades: List of clades to process
            
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
            
            # Extract beta values and substitution counts for each clade
            for clade in clades:
                site_dict[f'beta_{clade}'] = float(row[beta_idx_map[clade]])
                site_dict[f'subs_{clade}'] = float(row[subs_idx_map[clade]])
            
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
    def get_site_fields(clades: List[str]) -> List[str]:
        """Get list of site-specific fields produced by this method.
        
        Args:
            clades: List of clades to generate fields for
            
        Returns:
            List of field names
        """
        fields = []
        for clade in clades:
            fields.extend([
                f'beta_{clade}',
                f'subs_{clade}'
            ])
        return fields
