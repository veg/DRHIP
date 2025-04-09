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
from typing import Dict, List, Any

from ..utils import file_handlers as fh

# Define a lock for synchronizing writes to the output files
write_lock = threading.Lock()

def get_omega3(fit: dict) -> Dict[str, float]:
    """Get the omega3 value and proportion from a model fit.
    
    Args:
        fit: Dictionary containing model fit results
        
    Returns:
        Dictionary containing omega and proportion values from the rate distribution
    """
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    omega_idx = str(len(omegas) - 1)
    return omegas[omega_idx]


def process_gene(gene: str, results_path: str, output_dir: str, site_mappings_dir: str, 
                summary_fieldnames: List[str], sites_fieldnames: List[str]) -> None:
    """Process HyPhy results for a single gene.
    
    Args:
        gene: Name of the gene to process
        results_path: Path to the directory containing HyPhy results
        output_dir: Directory to write output files
        site_mappings_dir: Directory containing site mapping files
        summary_fieldnames: List of field names for summary CSV
        sites_fieldnames: List of field names for sites CSV
        
    The function processes multiple HyPhy analysis results:
    - RELAX: Tests for relaxation of selection
    - BUSTED: Branch-site unrestricted statistical test for episodic diversification
    - CFEL: Contrast-FEL for comparing selection between clades
    - FEL: Fixed effects likelihood test for selection
    - MEME: Mixed effects model of evolution
    - PRIME: Property Informed Models of Evolution
    """
    ####### Set up output files #######
    outfile_summary = os.path.join(output_dir, f"{gene}_summary.csv")
    outfile_sites = os.path.join(output_dir, f"{gene}_sites.csv")

    with write_lock:
        with open(outfile_summary, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=summary_fieldnames)
            writer.writeheader()

        with open(outfile_sites, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=sites_fieldnames)
            writer.writeheader()
    ####### End set up output files #######

    # Load RELAX JSON
    relax_path = os.path.join(results_path, "concat", "RELAX", f"{gene}.RELAX.json")
    relax = fh.load_json(relax_path)

    # Load BUSTED JSON
    busted_path = os.path.join(results_path, "concat", "BUSTED", f"{gene}.BUSTED.json")
    busted = fh.load_json(busted_path)

    if not busted or not relax:
        print(f"Missing required analysis files for {gene}")
        return

    # Load CFEL JSON
    cfel_path = os.path.join(results_path, "concat", "contrastFEL", f"{gene}.CFEL.json")
    cfel = fh.load_json(cfel_path)

    if not cfel:
        print(f"Missing CFEL analysis for {gene}")
        return

    ####### Preprocess CFEL data #######
    by_type: Dict[str, List[str]] = {}

    # Get the tag for each branch
    tested = cfel['tested']['0']

    # Group branch names by tag
    for branch, tag in tested.items():
        if tag not in by_type:
            by_type[tag] = []
        by_type[tag].append(branch)

    # Extract rows and headers from the CFEL data
    data_rows = cfel["MLE"]["content"]["0"]
    headers = cfel["MLE"]["headers"]

    # Build column lookup maps
    header_map = {short_label: i for i, (short_label, _) in enumerate(headers)}
    beta_idx_map = {c: header_map[f"beta ({c})"] for c in fh.all_clades}
    subs_idx_map = {c: header_map[f"subs ({c})"] for c in fh.all_clades}
    ####### END Preprocess CFEL data #######

    for clade in fh.all_clades:
        print(f"Processing {gene} in clade {clade}...")
        gene_summary_dict = {'gene': gene, 'clade': clade}

        invariant = [0, 0]  # [nt_conserved, aa_conserved]
        site_recorder: Dict[int, Dict[str, Any]] = {}

        ####### Process RELAX results #######
        # Get test results
        test_results = relax['test results']
        gene_summary_dict['relax_pvalue'] = test_results['p-value']
        gene_summary_dict['relax_k'] = test_results['relaxation or intensification parameter']
        gene_summary_dict['relax_lrt'] = test_results['LRT']
        gene_summary_dict['relax_k_significant'] = test_results['p-value'] <= 0.05

        ####### Process BUSTED results #######
        # Get test results
        busted_results = busted['test results']
        gene_summary_dict['busted_pvalue'] = busted_results['p-value']
        gene_summary_dict['busted_lrt'] = busted_results['LRT']
        gene_summary_dict['busted_evidence'] = busted_results['p-value'] <= 0.05

        # Get omega3 distribution
        omega3 = get_omega3(busted)
        gene_summary_dict['busted_omega3'] = omega3['omega']
        gene_summary_dict['busted_omega3_weight'] = omega3['weight']

        # Get branch lengths
        branch_lengths = busted['branch attributes']['0']
        total_branch_length = sum(float(branch['length']) for branch in branch_lengths.values())
        gene_summary_dict['total_branch_length'] = total_branch_length

        # Calculate tree length ratio
        if 'reference' in busted['branch attributes']:
            ref_branch_lengths = busted['branch attributes']['reference']
            ref_total_length = sum(float(branch['length']) for branch in ref_branch_lengths.values())
            gene_summary_dict['tree_length_ratio'] = total_branch_length / ref_total_length
        else:
            gene_summary_dict['tree_length_ratio'] = 1.0

        ####### Process site-specific analyses #######
        # Process each site in the data rows
        for row in data_rows:
            site = int(row[0])  # First column is the site index
            site_dict = {'gene': gene, 'site': site, 'clade': clade}

            # Extract beta values and substitution counts for each clade
            for c in fh.all_clades:
                site_dict[f'beta_{c}'] = float(row[beta_idx_map[c]])
                site_dict[f'subs_{c}'] = float(row[subs_idx_map[c]])

            # Check for conservation
            if all(site_dict[f'subs_{c}'] == 0 for c in fh.all_clades):
                invariant[0] += 1  # nucleotide conserved
                if all(site_dict[f'beta_{c}'] == 1 for c in fh.all_clades):
                    invariant[1] += 1  # amino acid conserved

            # Store site data for later use
            site_recorder[site] = site_dict

        # Record conservation statistics
        gene_summary_dict['nt_conserved'] = invariant[0]
        gene_summary_dict['aa_conserved'] = invariant[1]

        # Write results to files
        with write_lock:
            with open(outfile_summary, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=summary_fieldnames)
                writer.writerow(gene_summary_dict)

            with open(outfile_sites, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=sites_fieldnames)
                for site_dict in site_recorder.values():
                    writer.writerow(site_dict)
