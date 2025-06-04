"""
Core analysis logic for processing HyPhy results for individual genes.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import os
import csv
import threading
from typing import Dict, Any, Optional

from ..config import SUMMARY_FIELDNAMES, SITES_FIELDNAMES, DEFAULT_COMPARISON_GROUPS
from ..utils import file_handlers as fh
from ..utils.result_helpers import merge_method_data, ensure_ordered_fields, detect_comparison_groups
from ..methods import HyPhyMethodRegistry

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
    
    # Expected fields from config
    expected_summary_fields = set(SUMMARY_FIELDNAMES)
    expected_site_fields = set(SITES_FIELDNAMES)
    
    # Fields get added to the output as we process the methods
    output_summary_fields = {'gene'}  # Gene is always provided
    output_site_fields = {'gene', 'site'}  # These are always provided
    
    ####### Set up output files #######
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    outfile_summary = os.path.join(output_dir, f"{gene}_summary.csv")
    outfile_sites = os.path.join(output_dir, f"{gene}_sites.csv")
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
    
    # Check if we have at least one method with results
    if not method_results:
        print(f"No method results found for {gene}, skipping")
        return
        
    # Log which methods are missing but we're continuing anyway
    if missing_methods:
        print(f"Processing {gene} with {len(method_results)} methods. Missing: {', '.join(missing_methods)}")        
    
    # Detect comparison groups from available methods
    comparison_groups, groups_by_method = detect_comparison_groups(
        method_results, 
        default_groups=DEFAULT_COMPARISON_GROUPS,
        methods_to_check=['CFEL', 'RELAX']  # Can be extended with other methods in the future
    )
        
    # Make detected comparison groups available to methods
    for method in methods:
        if hasattr(method, 'set_comparison_groups'):
            method.set_comparison_groups(comparison_groups)
    

    # Initialize summary dictionary with gene name only
    gene_summary_dict = {'gene': gene}
    
    # Dictionary to track site-specific data
    site_recorder: Dict[int, Dict[str, Any]] = {}
    
    # Track which methods provide which fields
    field_providers = {}
    site_field_providers = {}

    # Process each method's results
    for method in methods:
        if method.name in method_results:
            # Process method-specific summary results
            summary_data = method.process_results(method_results[method.name])
            
            # Track which summary fields are provided by this method
            output_summary_fields.update(summary_data.keys())
            
            # Track which methods provide which fields
            for field in summary_data.keys():
                if field not in field_providers:
                    field_providers[field] = []
                field_providers[field].append(method.name)
            
            # Use helper to merge summary data
            merge_method_data(
                target_dict=gene_summary_dict,
                method_data=summary_data,
                method_name=method.name,
                providers=field_providers,
                context={'gene': gene}
            )

            # Process method-specific site data if available
            if hasattr(method, 'process_site_data'):
                site_data = method.process_site_data(method_results[method.name])
                for site, data in site_data.items():
                    if site not in site_recorder:
                        # Initialize site data with only gene and site
                        site_recorder[site] = {
                            'gene': gene, 
                            'site': site
                        }
                        
                    # Track fields that will be in the output
                    output_site_fields.update(data.keys())
                    
                    # Use helper to merge site data
                    merge_method_data(
                        target_dict=site_recorder[site],
                        method_data=data,
                        method_name=method.name,
                        providers=site_field_providers,
                        context={'gene': gene, 'site': site}
                    )

    # Validate that all expected summary fields are present
    if not expected_summary_fields.issubset(output_summary_fields):
        missing_summary_fields = expected_summary_fields - output_summary_fields
        print(f"Warning: Missing summary fields for {gene}: {missing_summary_fields}")
        
    # Validate that all expected site fields are present
    if not expected_site_fields.issubset(output_site_fields):
        missing_site_fields = expected_site_fields - output_site_fields
        print(f"Warning: Missing site fields for {gene}: {missing_site_fields}")

    # Check if files exist and warn about overwriting
    if os.path.exists(outfile_summary):
        print(f"Warning: Overwriting existing file {outfile_summary}")
    if os.path.exists(outfile_sites):
        print(f"Warning: Overwriting existing file {outfile_sites}")
    
    # Use helper to ensure fields are ordered properly
    ordered_summary_fields = ensure_ordered_fields(output_summary_fields, ['gene'])
    ordered_site_fields = ensure_ordered_fields(output_site_fields, ['gene', 'site'])
    
    # Write results to files
    with write_lock:
        # Write summary data
        with open(outfile_summary, 'w', newline='') as csvfile:  # 'w' mode to overwrite
            writer = csv.DictWriter(csvfile, fieldnames=ordered_summary_fields)
            writer.writeheader()  # Always write header
            writer.writerow(gene_summary_dict)

        # Write site data
        with open(outfile_sites, 'w', newline='') as csvfile:  # 'w' mode to overwrite
            writer = csv.DictWriter(csvfile, fieldnames=ordered_site_fields, extrasaction='ignore')
            writer.writeheader()  # Always write header
            for site_dict in site_recorder.values():
                writer.writerow(site_dict)
