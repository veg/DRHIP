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
