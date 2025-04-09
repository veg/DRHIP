"""
Registry for HyPhy analysis methods.
"""

from typing import Dict, List

from .base import HyPhyMethod
from .busted import BustedMethod
from .relax import RelaxMethod
from .cfel import CfelMethod
from .fel import FelMethod
from .meme import MemeMethod
from .prime import PrimeMethod

class HyPhyMethodRegistry:
    """Registry for HyPhy analysis methods."""
    
    def __init__(self):
        """Initialize an empty registry."""
        self._methods: Dict[str, HyPhyMethod] = {}
        
        # Register built-in methods
        self.register(BustedMethod())
        self.register(RelaxMethod())
        self.register(CfelMethod())
        self.register(FelMethod())
        self.register(MemeMethod())
        self.register(PrimeMethod())
    
    def register(self, method: HyPhyMethod) -> None:
        """Register a new method.
        
        Args:
            method: HyPhyMethod instance to register
        """
        self._methods[method.name] = method
    
    def get_method(self, name: str) -> HyPhyMethod:
        """Get a registered method by name.
        
        Args:
            name: Name of the method to retrieve
            
        Returns:
            The requested HyPhyMethod instance
            
        Raises:
            KeyError: If method is not registered
        """
        return self._methods[name]
    
    def get_all_methods(self) -> List[HyPhyMethod]:
        """Get all registered methods.
        
        Returns:
            List of all registered HyPhyMethod instances
        """
        return list(self._methods.values())
    
    def get_all_summary_fields(self) -> List[str]:
        """Get all summary fields from all registered methods.
        
        Returns:
            List of all summary field names
        """
        fields = set()
        for method in self._methods.values():
            if hasattr(method, 'get_summary_fields'):
                fields.update(method.get_summary_fields())
        return sorted(list(fields))
    
    def get_all_site_fields(self, clades: List[str]) -> List[str]:
        """Get all site-specific fields from all registered methods.
        
        Args:
            clades: List of clades to get fields for
            
        Returns:
            List of all site field names
        """
        fields = set()
        for method in self._methods.values():
            if hasattr(method, 'get_site_fields'):
                fields.update(method.get_site_fields(clades))
        return sorted(list(fields))
