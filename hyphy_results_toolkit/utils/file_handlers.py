"""
File handling utilities for HyPhy results.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import os
import json


def get_genes(results_path: str) -> list:
    """Get list of genes from the results directory.
    
    Args:
        results_path: Path to the directory containing HyPhy results
        
    Returns:
        List of gene names found in the results directory
    """
    concat_dir = os.path.join(results_path, "concat")
    if not os.path.exists(concat_dir):
        raise FileNotFoundError(f"Expected to find directory 'concat' in {results_path}")
    
    busted_dir = os.path.join(concat_dir, "BUSTED")
    if not os.path.exists(busted_dir):
        raise FileNotFoundError(f"Expected to find directory 'BUSTED' in {concat_dir}")
    
    genes = []
    for file in os.listdir(busted_dir):
        if file.endswith(".BUSTED.json"):
            genes.append(file.split(".BUSTED.json")[0])
    
    return sorted(genes)


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
