"""
Base class and interfaces for HyPhy analysis methods.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, List

class HyPhyMethod(ABC):
    """Base class for HyPhy analysis methods."""
    
    def __init__(self, name: str, file_suffix: str):
        """Initialize a HyPhy method.
        
        Args:
            name: Name of the method (e.g., 'BUSTED', 'RELAX')
            file_suffix: Suffix for result files (e.g., 'BUSTED.json')
        """
        self.name = name
        self.file_suffix = file_suffix
        self._comparison_groups = None
        
    def set_comparison_groups(self, comparison_groups):
        """Set the comparison groups to use for analysis.
        
        Args:
            comparison_groups: List of comparison group names
        """
        self._comparison_groups = comparison_groups
        return self
    
    def get_file_path(self, results_path: str, gene: str) -> str:
        """Get the path to results file for this method.
        
        Args:
            results_path: Base path to results directory
            gene: Name of the gene
            
        Returns:
            Full path to the results file
        """
        import os
        from ..config import METHOD_PATHS
        method_dir = METHOD_PATHS.get(self.name, self.name)
        return os.path.join(results_path, method_dir, f"{gene}.{self.file_suffix}")
    
    @abstractmethod
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process method-specific results.
        
        Args:
            results: Raw results dictionary from JSON file
            
        Returns:
            Processed results dictionary with standardized keys
        """
        pass
    
    def process_site_data(self, site_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific data.
        
        Args:
            site_data: Site-specific data from results
            
        Returns:
            Processed site data dictionary with fields that are consistent across all sites
        """
        return {}
        
    def process_comparison_site_data(self, site_data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data.
        
        Args:
            site_data: Site-specific data from results
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data.
            The inner dictionaries map comparison group names to field values.
            Format: {site_id: {group_name: {field: value}}}
        """
        return {}
        
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []
        
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return []
        
    @staticmethod
    def get_comparison_group_fields() -> List[str]:
        """Get list of fields that are specific to comparison groups.
        
        These fields will be written to a separate CSV file with comparison_group column.
        """
        return []
    
    def has_mle_content(self, results: Dict[str, Any]) -> bool:
        """Check if results have MLE content structure.
        
        Args:
            results: Raw results dictionary
            
        Returns:
            True if MLE content structure is present, False otherwise
        """
        return ('MLE' in results and 
                'content' in results.get('MLE', {}) and 
                '0' in results.get('MLE', {}).get('content', {}))
    
    def has_mle_headers(self, results: Dict[str, Any]) -> bool:
        """Check if results have MLE headers structure.
        
        Args:
            results: Raw results dictionary
            
        Returns:
            True if MLE headers structure is present, False otherwise
        """
        return ('MLE' in results and 
                'headers' in results.get('MLE', {}))
    
    def get_header_indices(self, results: Dict[str, Any]) -> Dict[str, int]:
        """Get mapping from header names to column indices.
        
        Args:
            results: Raw results dictionary
            
        Returns:
            Dictionary mapping header names to column indices
        """
        header_indices = {}
        
        if self.has_mle_headers(results):
            headers = results['MLE']['headers']
            for i, (header_name, _) in enumerate(headers):
                header_indices[header_name] = i
        
        return header_indices
    
    def get_column_index(self, header_indices: Dict[str, int], 
                         column_name: str, default_index: int) -> int:
        """Get index for a specific column with fallback to default.
        
        Args:
            header_indices: Dictionary mapping header names to column indices
            column_name: Name of the column to find
            default_index: Default index to use if column name not found
            
        Returns:
            Column index
        """
        return header_indices.get(column_name, default_index)

    def process_site_mle_data(self, results: Dict[str, Any], column_names: Dict[str, str], 
                           process_row_fn) -> Dict[int, Dict[str, Any]]:
        """Process site-specific data from MLE content.
        
        This helper handles common site data processing patterns:
        - Checking for valid MLE content and headers
        - Extracting column indices
        - Processing each site with error handling
        - Returning a dictionary of site-specific results
        
        Args:
            results: Raw results dictionary from JSON file
            column_names: Dictionary mapping logical column names to actual column names in the results
                         e.g. {'alpha': 'alpha', 'beta': 'beta', 'p-value': 'p-value'}
            process_row_fn: Function that takes (site_idx, row, column_indices) and returns site data dict
            
        Returns:
            Dictionary mapping site indices to site-specific data
        """
        site_results = {}
        
        # Check if required data is available
        if not self.has_mle_content(results) or not self.has_mle_headers(results):
            return site_results
            
        # Get header indices
        header_indices = self.get_header_indices(results)
        
        # Get indices for the values we need
        column_indices = {}
        for logical_name, actual_name in column_names.items():
            # Default to -1 to indicate column not found
            column_indices[logical_name] = self.get_column_index(header_indices, actual_name, -1)
        
        # Process each site (row index + 1 is the site number)
        for site_idx, row in enumerate(results['MLE']['content']['0'], 1):
            try:
                # Process the row using the provided function
                site_data = process_row_fn(site_idx, row, column_indices)
                if site_data:
                    site_results[site_idx] = site_data
            except Exception as e:
                print(f"Error processing site {site_idx}: {e}")
                # Skip this site on error
                continue
        
        return site_results
    
    def extract_common_fields(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract common fields from HyPhy results.
        
        This method extracts fields that are common across multiple HyPhy methods:
        - N: Number of sequences
        - T: Total branch length
        - sites: Number of sites
        
        Args:
            results: Raw results dictionary from JSON file
            
        Returns:
            Dictionary with extracted common fields
        """
        common_fields = {}
        
        # Extract number of sites
        if 'input' in results:
            # Try different keys that might contain the number of sites
            if 'number of sites' in results['input']:
                common_fields['sites'] = results['input']['number of sites']
            elif 'sites' in results['input']:
                common_fields['sites'] = results['input']['sites']
        
        # If sites not found in input, try to get from MLE content length
        if 'sites' not in common_fields and self.has_mle_content(results):
            common_fields['sites'] = len(results['MLE']['content']['0'])
        

        # Get tested branches
        if 'tested' in results and '0' in results['tested']:
            tested_branches = {bn: val for bn,val in results['tested']['0'].items() if val == "test"}

            # Extract number of tested branches
            common_fields['N'] = len(tested_branches)
        
        # Extract total branch length
        if 'branch attributes' in results and '0' in results['branch attributes']:
            try:
                branch_lengths = []
                branch_data = results['branch attributes']['0']

                # keep only tested branches
                if tested_branches:
                    branch_data = {bn: val for bn,val in branch_data.items() if bn in tested_branches}
                
                # Handle different formats of branch length data
                for branch in branch_data.values():
                    # Format 1: Direct 'length' key
                    if isinstance(branch, dict) and 'length' in branch:
                        branch_lengths.append(float(branch['length']))
                    # Format 2: MG94xREV key
                    elif isinstance(branch, dict) and 'MG94xREV with separate rates for branch sets' in branch:
                        branch_lengths.append(float(branch['MG94xREV with separate rates for branch sets']))
                    # Format 3: Direct float value
                    elif isinstance(branch, (int, float, str)):
                        try:
                            branch_lengths.append(float(branch))
                        except (ValueError, TypeError):
                            pass
                
                if branch_lengths:
                    common_fields['T'] = sum(branch_lengths)
            except (ValueError, TypeError, KeyError):
                pass
        
        return common_fields
