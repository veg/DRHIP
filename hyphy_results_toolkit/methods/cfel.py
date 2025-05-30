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
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific CFEL data.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        site_results = {}
        
        # Get the tag for each branch
        tested = results['tested']['0']
        by_type = {}
        
        # Group branch names by tag
        for branch, tag in tested.items():
            if tag not in by_type:
                by_type[tag] = []
            by_type[tag].append(branch)
        
        # Build column lookup maps
        headers = results["MLE"]["headers"]
        self._build_column_maps(headers, list(by_type.keys()))
        
        # Process site data
        data_rows = results["MLE"]["content"]["0"]
        p_value_idx = self._header_map.get("p-value", 5)  # Default to index 5 if not found
        
        for row in data_rows:
            site = int(row[0])  # First column is the site index
            p_value = float(row[p_value_idx])
            
            # Create a formatted marker string for the cfel_marker field
            if p_value <= 0.05:
                # Format similar to End2End-DENV
                comparisons = []
                comparisons.append(f"overall: {p_value:.3f}")
                
                # Add pairwise comparisons
                for i, group1 in enumerate(by_type.keys()):
                    for group2 in list(by_type.keys())[i+1:]:
                        comparison_key = f"{group1} vs {group2}"
                        comparisons.append(f"{comparison_key}: {p_value:.3f}")
                        
                marker = ", ".join(comparisons)
            else:
                marker = "-"
            
            site_results[site] = {
                'cfel_marker': marker
            }
        
        return site_results
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []  # No summary fields needed
    
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'cfel_marker'
        ]