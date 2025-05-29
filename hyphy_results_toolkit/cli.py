"""
Command-line interface for hyphy-results-toolkit.

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import argparse
import os
import concurrent.futures
import traceback

from .parsers import process_gene
from .utils import file_handlers as fh


def main():
    """Main entry point for the CLI."""
    arguments = argparse.ArgumentParser(
        description='Summarize HyPhy analysis results for many genes into two csv files'
    )
    arguments.add_argument(
        '-i', '--input',
        help='Path to hyphy results directory (CAPHEINE workflow format)',
        required=False,
        default=os.getcwd(),
        type=str
    )
    arguments.add_argument(
        '-o', '--output',
        help='Path to output directory (defaults to current directory)',
        required=False,
        default=os.getcwd(),
        type=str
    )
    settings = arguments.parse_args()

    # Handle both absolute and relative paths for hyphy_results
    if os.path.isabs(settings.input):
        results_path = settings.input
    else:
        results_path = os.path.join(os.getcwd(), settings.input)
    
    # Use output directory from arguments
    output_dir = settings.output
    
    # Check if site_mappings directory exists
    site_mappings_dir = os.path.join(results_path, "site_mappings")
    if not os.path.exists(site_mappings_dir):
        print(f"Warning: site_mappings directory not found at {site_mappings_dir}")
        site_mappings_dir = None

    genes = fh.get_genes(results_path)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                process_gene.process_gene,
                gene,
                results_path,
                output_dir,
                site_mappings_dir
            ) for gene in genes
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                tb = traceback.format_exc()
                print(f'Generated an exception: {exc}\nTraceback: {tb}')


if __name__ == "__main__":
    main()
