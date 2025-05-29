"""
File handling utilities for HyPhy results.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import os
import json
from typing import List

from ..methods.registry import HyPhyMethodRegistry
from ..config import METHOD_PATHS

# Default comparison groups 
# These represent different sequence groups being compared in selection analyses
comparison_groups = ['foreground', 'background']

def get_genes(results_path: str) -> List[str]:
    """Get list of genes from the results directory.
    
    Args:
        results_path: Path to the directory containing HyPhy results
        
    Returns:
        List of gene names found in the results directory across all methods
    """
    # Get methods from registry
    registry = HyPhyMethodRegistry()
    methods = registry.get_all_methods()
    
    # Set to store unique gene names
    all_genes = set()
    
    # Collect genes from all method directories
    for method in methods:
        # Get method name and file suffix
        method_name = method.name
        file_suffix = method.file_suffix
        
        # Get method directory from config
        method_dir_name = METHOD_PATHS.get(method_name, method_name)
        method_dir = os.path.join(results_path, method_dir_name)
        
        if os.path.exists(method_dir) and os.path.isdir(method_dir):
            for file in os.listdir(method_dir):
                # Extract gene name from file name based on method suffix
                if file.endswith(f".{file_suffix}"):
                    gene = file.split(f".{file_suffix}")[0]
                    all_genes.add(gene)
    
    if not all_genes:
        raise FileNotFoundError(f"No HyPhy analysis results found in {results_path}")
    
    return sorted(list(all_genes))


def load_json(filepath: str) -> dict:
    """Load and parse a JSON file.
    
    Args:
        filepath: Path to JSON file
        
    Returns:
        Parsed JSON content as dictionary
    """
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    except json.JSONDecodeError:
        print(f"Error decoding JSON from file: {filepath}")
        return None
