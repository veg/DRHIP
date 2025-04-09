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


def process_gene(gene: str, results_path: str, output_dir: str, site_mappings_dir: str) -> None:
    """Process HyPhy results for a single gene.
    
    Args:
        gene: Name of the gene to process
        results_path: Path to the directory containing HyPhy results
        output_dir: Directory to write output files
        site_mappings_dir: Directory containing site mapping files
        
    The function processes multiple HyPhy analysis results using registered methods.
    Each method processes its own results and contributes to summary and site-specific data.
    """
    # Initialize method registry
    registry = HyPhyMethodRegistry()
    
    # Get all methods
    methods = registry.get_all_methods()
    
    # Get complete field lists
    summary_fields = SUMMARY_FIELDNAMES + registry.get_all_summary_fields()
    site_fields = SITES_FIELDNAMES + registry.get_all_site_fields(fh.all_clades)
    
    ####### Set up output files #######
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
    for method in methods:
        file_path = method.get_file_path(results_path, gene)
        result = fh.load_json(file_path)
        if result:
            method_results[method.name] = result
        else:
            print(f"Missing {method.name} analysis for {gene}")
            if method.name in ['BUSTED', 'RELAX', 'CFEL']:  # Required methods
                return

    for clade in fh.all_clades:
        print(f"Processing {gene} in clade {clade}...")
        gene_summary_dict = {'gene': gene, 'clade': clade}
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
                            site_recorder[site] = {'gene': gene, 'site': site, 'clade': clade}
                        site_recorder[site].update(data)

        # Write results to files
        with write_lock:
            with open(outfile_summary, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=summary_fields)
                writer.writerow(gene_summary_dict)

            with open(outfile_sites, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=site_fields)
                for site_dict in site_recorder.values():
                    writer.writerow(site_dict)
