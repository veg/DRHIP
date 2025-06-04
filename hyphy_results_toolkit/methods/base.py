"""
Base class and interfaces for HyPhy analysis methods.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any

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
            Processed site data dictionary
        """
        return {}
    
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

    # TODO: some of this looks like it might have been ai hallucinations lol
    # but it also looks to have produced the right values so far
    # i suspect maybe its needlessly complicated is all, for checking a bit of nonsense
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
        
        # Extract number of sequences
        if 'input' in results:
            # Try different keys that might contain the number of sequences
            if 'number of sequences' in results['input']:
                common_fields['N'] = results['input']['number of sequences']
            elif 'sequences' in results['input']:
                common_fields['N'] = len(results['input']['sequences'])
        
        # If N not found in input, try to get from tested sequences
        if 'N' not in common_fields and 'tested' in results and 'sequences' in results['tested']:
            common_fields['N'] = len(results['tested']['sequences'])
        
        # Extract total branch length
        if 'branch attributes' in results and '0' in results['branch attributes']:
            try:
                branch_lengths = []
                branch_data = results['branch attributes']['0']
                
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
