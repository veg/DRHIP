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
from typing import Dict, Any, Optional, Tuple

from ..config import SUMMARY_FIELDNAMES, SITES_FIELDNAMES, DEFAULT_COMPARISON_GROUPS
from ..utils import file_handlers as fh
from ..utils.result_helpers import (
    merge_method_data,
    ensure_ordered_fields,
    detect_comparison_groups,
    collect_method_fields,
    validate_fields,
)
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
    # Track expected field names for validation
    expected_summary_fields = set(SUMMARY_FIELDNAMES)
    expected_site_fields = set(SITES_FIELDNAMES)

    # Track which fields are actually in the output
    output_summary_fields = set(["gene"])  # Always include gene
    output_site_fields = set(["gene", "site"])  # Always include gene and site
    output_comparison_site_fields = set(["gene", "site", "comparison_group"])
    output_comparison_summary_fields = set(
        ["gene", "comparison_group"]
    )  # For non-site-specific comparison data

    ####### Set up output files #######
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    outfile_summary = os.path.join(output_dir, f"{gene}_summary.csv")
    outfile_sites = os.path.join(output_dir, f"{gene}_sites.csv")
    outfile_comparison_site = os.path.join(output_dir, f"{gene}_comparison_site.csv")
    outfile_comparison_summary = os.path.join(
        output_dir, f"{gene}_comparison_summary.csv"
    )
    ####### End set up output files #######

    # Load and validate method results
    method_results: Dict[str, Optional[Dict]] = {}
    missing_methods = []

    for method in methods:
        file_path = method.get_file_path(results_path, gene)
        result = fh.load_json(file_path)
        if result:
            # Per-method JSON validation
            missing_fields = []
            if hasattr(method, "validate_input_json"):
                try:
                    missing_fields = method.validate_input_json(result) or []
                except Exception as e:
                    print(f"Error during validation for {method.name} on {gene}: {e}")
                    missing_fields = ["<validation error>"]

            if missing_fields:
                print(
                    f"WARNING: {method.name} for {gene} is missing required fields: {missing_fields}. Final output may be incomplete. Check that you're using the correct version of HyPhy."
                )

            method_results[method.name] = result
        else:
            missing_methods.append(method.name)

    # Check if we have at least one method with results
    if not method_results:
        print(f"No method results found for {gene}, skipping")
        return

    # Log which methods are missing but we're continuing anyway
    if missing_methods:
        print(
            f"Processing {gene} with {len(method_results)} methods. Missing: {', '.join(missing_methods)}"
        )

    # Detect comparison groups from available methods
    comparison_groups, groups_by_method = detect_comparison_groups(
        method_results,
        default_groups=DEFAULT_COMPARISON_GROUPS,
        methods_to_check=[
            "CFEL",
            "RELAX",
        ],  # Can be extended with other methods in the future
        gene_name=gene,
    )

    # Make detected comparison groups available to methods
    for method in methods:
        if hasattr(method, "set_comparison_groups"):
            method.set_comparison_groups(comparison_groups)

    # Initialize summary dictionary with gene name only
    gene_summary_dict = {"gene": gene}

    # Dictionary to track site-specific data
    site_recorder: Dict[int, Dict[str, Any]] = {}

    # Dictionary to track site-specific comparison group data
    comparison_site_recorder: Dict[Tuple[int, str], Dict[str, Any]] = {}

    # Dictionary to track non-site-specific comparison group data
    comparison_summary_recorder: Dict[str, Dict[str, Any]] = {}

    # Track which methods provide which fields
    field_providers = {}
    site_field_providers = {}
    comparison_site_field_providers = {}
    comparison_summary_field_providers = {}

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
                context={"gene": gene},
            )

            # Process method-specific site data if available
            if hasattr(method, "process_site_data"):
                site_data = method.process_site_data(method_results[method.name])
                for site, data in site_data.items():
                    if site not in site_recorder:
                        # Initialize site data with only gene and site
                        site_recorder[site] = {"gene": gene, "site": site}

                    # Track fields that will be in the output
                    output_site_fields.update(data.keys())

                    # Use helper to merge site data
                    merge_method_data(
                        target_dict=site_recorder[site],
                        method_data=data,
                        method_name=method.name,
                        providers=site_field_providers,
                        context={"gene": gene, "site": site},
                    )

            # Process comparison group-specific site data if available
            if comparison_groups and hasattr(method, "process_comparison_site_data"):
                comparison_site_data = method.process_comparison_site_data(
                    method_results[method.name]
                )

                # Process each site's comparison group data
                for site, site_comp_data in comparison_site_data.items():
                    for group, group_data in site_comp_data.items():
                        # Create a unique key for each site-group combination
                        comp_key = (site, group)
                        if comp_key not in comparison_site_recorder:
                            comparison_site_recorder[comp_key] = {
                                "gene": gene,
                                "site": site,
                                "comparison_group": group,
                            }

                        # Track fields that will be in the output
                        output_comparison_site_fields.update(group_data.keys())

                        # Use helper to merge comparison group data
                        merge_method_data(
                            target_dict=comparison_site_recorder[comp_key],
                            method_data=group_data,
                            method_name=method.name,
                            providers=comparison_site_field_providers,
                            context={
                                "gene": gene,
                                "site": site,
                                "comparison_group": group,
                            },
                        )

            # Process comparison group-specific non-site data if available
            if comparison_groups and hasattr(method, "process_comparison_data"):
                comparison_data = method.process_comparison_data(
                    method_results[method.name]
                )

                # Process each group's non-site-specific data
                for group, group_data in comparison_data.items():
                    # Use group name as key
                    if group not in comparison_summary_recorder:
                        comparison_summary_recorder[group] = {
                            "gene": gene,
                            "comparison_group": group,
                        }

                    # Track fields that will be in the output
                    output_comparison_summary_fields.update(group_data.keys())

                    # Use helper to merge comparison group data
                    merge_method_data(
                        target_dict=comparison_summary_recorder[group],
                        method_data=group_data,
                        method_name=method.name,
                        providers=comparison_summary_field_providers,
                        context={"gene": gene, "comparison_group": group},
                    )

    # Validate summary and site fields
    validate_fields(
        expected_summary_fields,
        output_summary_fields,
        context="summary",
        entity_name=gene,
    )
    validate_fields(
        expected_site_fields, output_site_fields, context="site", entity_name=gene
    )

    # Only validate comparison fields if we have comparison groups
    if comparison_groups:
        # Collect and validate site-specific comparison fields from active methods
        active_comparison_site_fields = collect_method_fields(
            methods=methods,
            method_results=method_results,
            field_getter_name="get_comparison_group_site_fields",
            base_fields={"gene", "site", "comparison_group"},
        )
        validate_fields(
            active_comparison_site_fields,
            output_comparison_site_fields,
            context="comparison site",
            entity_name=gene,
        )

        # Collect and validate non-site-specific comparison fields from active methods
        active_comparison_summary_fields = collect_method_fields(
            methods=methods,
            method_results=method_results,
            field_getter_name="get_comparison_group_summary_fields",
            base_fields={"gene", "comparison_group"},
        )
        validate_fields(
            active_comparison_summary_fields,
            output_comparison_summary_fields,
            context="comparison summary",
            entity_name=gene,
        )

    # Check if files exist and warn about overwriting
    if os.path.exists(outfile_summary):
        print(f"Warning: Overwriting existing file {outfile_summary}")
    if os.path.exists(outfile_sites):
        print(f"Warning: Overwriting existing file {outfile_sites}")

    # Use helper to ensure fields are ordered properly
    ordered_summary_fields = ensure_ordered_fields(output_summary_fields, ["gene"])
    ordered_site_fields = ensure_ordered_fields(output_site_fields, ["gene", "site"])
    ordered_comparison_site_fields = ensure_ordered_fields(
        output_comparison_site_fields, ["gene", "site", "comparison_group"]
    )
    ordered_comparison_summary_fields = ensure_ordered_fields(
        output_comparison_summary_fields, ["gene", "comparison_group"]
    )

    # Write results to files
    with write_lock:
        # Write summary data
        with open(outfile_summary, "w", newline="") as csvfile:  # 'w' mode to overwrite
            writer = csv.DictWriter(csvfile, fieldnames=ordered_summary_fields)
            writer.writeheader()  # Always write header
            writer.writerow(gene_summary_dict)

        # Write site data
        with open(outfile_sites, "w", newline="") as csvfile:  # 'w' mode to overwrite
            writer = csv.DictWriter(
                csvfile, fieldnames=ordered_site_fields, extrasaction="ignore"
            )
            writer.writeheader()  # Always write header
            for site_dict in site_recorder.values():
                writer.writerow(site_dict)

        # Write site-specific comparison group data if we have any
        if comparison_site_recorder:
            with open(
                outfile_comparison_site, "w", newline=""
            ) as csvfile:  # 'w' mode to overwrite
                writer = csv.DictWriter(
                    csvfile,
                    fieldnames=ordered_comparison_site_fields,
                    extrasaction="ignore",
                )
                writer.writeheader()  # Always write header
                for comp_dict in comparison_site_recorder.values():
                    writer.writerow(comp_dict)

        # Write non-site-specific comparison group data if we have any
        if comparison_summary_recorder:
            with open(
                outfile_comparison_summary, "w", newline=""
            ) as csvfile:  # 'w' mode to overwrite
                writer = csv.DictWriter(
                    csvfile,
                    fieldnames=ordered_comparison_summary_fields,
                    extrasaction="ignore",
                )
                writer.writeheader()  # Always write header
                for comp_dict in comparison_summary_recorder.values():
                    writer.writerow(comp_dict)
