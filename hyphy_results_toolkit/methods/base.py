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
        method_dir = METHOD_PATHS.get(self.name, f"concat/{self.name}")
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
