"""
Core analysis logic for processing HyPhy results for individual genes.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import os
import threading
import csv
from typing import Dict, Any, Optional

from ..utils import file_handlers as fh
from ..methods import HyPhyMethodRegistry
from ..config import SUMMARY_FIELDNAMES, SITES_FIELDNAMES

# Define a lock for synchronizing writes to the output files
write_lock = threading.Lock()


def process_gene(gene: str, results_path: str, output_dir: str) -> None:
    """Process HyPhy results for a single gene.
    
    Args:
        gene: Name of the gene to process
        results_path: Path to the directory containing HyPhy results
        output_dir: Directory to write output files
        
    The function processes multiple HyPhy analysis results using registered methods.
    Each method processes its own results and contributes to summary and site-specific data.
    """
    # Initialize method registry
    registry = HyPhyMethodRegistry()
    
    # Get all methods
    methods = registry.get_all_methods()
    
    # Get complete field lists
    summary_fields = SUMMARY_FIELDNAMES + registry.get_all_summary_fields()
    site_fields = SITES_FIELDNAMES + registry.get_all_site_fields(fh.comparison_groups)
    
    # Dictionary to track all site fields that might be added dynamically
    all_site_fields = set(site_fields)
    
    ####### Set up output files #######
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    outfile_summary = os.path.join(output_dir, f"{gene}_summary.csv")
    outfile_sites = os.path.join(output_dir, f"{gene}_sites.csv")

    with write_lock:
        with open(outfile_summary, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=summary_fields)
            writer.writeheader()

        with open(outfile_sites, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=site_fields)
            writer.writeheader()
    ####### End set up output files #######

    # Load and validate method results
    method_results: Dict[str, Optional[Dict]] = {}
    missing_methods = []
    
    for method in methods:
        file_path = method.get_file_path(results_path, gene)
        result = fh.load_json(file_path)
        if result:
            method_results[method.name] = result
        else:
            missing_methods.append(method.name)
            print(f"Missing {method.name} analysis for {gene}")
    
    # Check if we have at least one method with results
    if not method_results:
        print(f"No method results found for {gene}, skipping")
        return
        
    # Log which methods are missing but we're continuing anyway
    if missing_methods:
        print(f"Processing {gene} with {len(method_results)} methods. Missing: {', '.join(missing_methods)}")        

    for group in fh.comparison_groups:
        print(f"Processing {gene} in comparison group {group}...")
        gene_summary_dict = {'gene': gene, 'comparison_group': group}
        site_recorder: Dict[int, Dict[str, Any]] = {}

        # Process each method's results
        for method in methods:
            if method.name in method_results:
                # Process method-specific summary results
                summary_data = method.process_results(method_results[method.name])
                gene_summary_dict.update(summary_data)

                # Process method-specific site data if available
                if hasattr(method, 'process_site_data'):
                    site_data = method.process_site_data(method_results[method.name])
                    for site, data in site_data.items():
                        if site not in site_recorder:
                            site_recorder[site] = {'gene': gene, 'site': site, 'comparison_group': group}
                        
                        # Track any new fields that weren't in the original site_fields list
                        for field in data.keys():
                            all_site_fields.add(field)
                            
                        site_recorder[site].update(data)

        # Convert set to list for CSV writer
        complete_site_fields = list(all_site_fields)
                        
        # Write results to files
        with write_lock:
            with open(outfile_summary, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=summary_fields)
                writer.writerow(gene_summary_dict)

            # Check if the sites file exists and has content
            header_exists = os.path.exists(outfile_sites) and os.path.getsize(outfile_sites) > 0
                
            with open(outfile_sites, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=complete_site_fields, extrasaction='ignore')
                if not header_exists:
                    writer.writeheader()
                for site_dict in site_recorder.values():
                    writer.writerow(site_dict)
