"""
Command-line interface for DRHIP (Data Reduction for HyPhy with Inference Processing).

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Hannah Verdonk (hannah.verdonk@temple.edu)
    Danielle Callan (dcallan@temple.edu)
"""

import argparse
import os
import concurrent.futures
import traceback
import csv
import tempfile

from .parsers import process_gene
from .utils import file_handlers as fh


def combine_csv_files(temp_dir: str, output_dir: str, file_suffix: str) -> None:
    """Combine all gene-specific CSV files into a single file with the superset of columns.
    
    Args:
        temp_dir: Directory containing gene-specific CSV files
        output_dir: Directory to write the combined file
        file_suffix: Suffix of the files to combine ('summary' or 'sites')
    """
    # Find all files with the given suffix
    files_to_combine = []
    for file in os.listdir(temp_dir):
        # Make sure we're only matching the exact suffix pattern
        # This prevents comparison_summary files from matching when looking for summary files
        if file.endswith(f"_{file_suffix}.csv") and not (
            file_suffix == "summary" and "_comparison_summary.csv" in file or
            file_suffix == "sites" and "_comparison_site.csv" in file
        ):
            files_to_combine.append(os.path.join(temp_dir, file))
    
    if not files_to_combine:
        print(f"No {file_suffix} files found to combine")
        return
    
    # Collect all unique fieldnames across files
    all_fieldnames = set()
    for file_path in files_to_combine:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            try:
                header = next(reader)
                all_fieldnames.update(header)
            except StopIteration:
                # Skip empty files
                continue
    
    # Ensure 'gene' is the first column for summary files
    # For sites files, ensure 'gene' and 'site' are the first two columns
    # For comparison_site files, ensure 'gene', 'site', and 'comparison_group' are the first three columns
    # For comparison_summary files, ensure 'gene' and 'comparison_group' are the first two columns
    ordered_fieldnames = []
    if file_suffix == 'summary':
        ordered_fieldnames = ['gene']
        for field in all_fieldnames:
            if field != 'gene':
                ordered_fieldnames.append(field)
    elif file_suffix == 'sites':
        ordered_fieldnames = ['gene', 'site']
        for field in all_fieldnames:
            if field not in ['gene', 'site']:
                ordered_fieldnames.append(field)
    elif file_suffix == 'comparison_site':
        ordered_fieldnames = ['gene', 'site', 'comparison_group']
        for field in all_fieldnames:
            if field not in ['gene', 'site', 'comparison_group']:
                ordered_fieldnames.append(field)
    elif file_suffix == 'comparison_summary':
        ordered_fieldnames = ['gene', 'comparison_group']
        for field in all_fieldnames:
            if field not in ['gene', 'comparison_group']:
                ordered_fieldnames.append(field)
    
    # Create the combined output file
    output_file = os.path.join(output_dir, f"combined_{file_suffix}.csv")
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=ordered_fieldnames)
        writer.writeheader()
        
        # Read each input file and write rows to the combined file
        for file_path in files_to_combine:
            with open(file_path, 'r', newline='') as infile:
                reader = csv.DictReader(infile)
                for row in reader:
                    # Fill missing fields with 'NA'
                    for field in ordered_fieldnames:
                        if field not in row:
                            row[field] = 'NA'
                    writer.writerow(row)
    
    print(f"Created combined {file_suffix} file: {output_file}")


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
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Create a temporary directory for gene-specific files
    with tempfile.TemporaryDirectory() as temp_dir:
        genes = fh.get_genes(results_path)
        print(f"Processing {len(genes)} genes...")

        # Process genes in parallel, storing results in temp directory
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(
                    process_gene.process_gene,
                    gene,
                    results_path,
                    temp_dir
                ) for gene in genes
            ]
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    tb = traceback.format_exc()
                    print(f'Generated an exception: {exc}\nTraceback: {tb}')
        
        # Combine gene-specific files into unified files
        print("Combining gene-specific results into unified files...")
        combine_csv_files(temp_dir, output_dir, "summary")
        combine_csv_files(temp_dir, output_dir, "sites")
        combine_csv_files(temp_dir, output_dir, "comparison_site")
        combine_csv_files(temp_dir, output_dir, "comparison_summary")


if __name__ == "__main__":
    main()
