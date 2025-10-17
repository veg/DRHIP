"""
Helper functions for processing HyPhy results.

This module contains utility functions for merging and processing results
from different HyPhy methods.
"""

from typing import Any, Dict, List, Optional, Set, Tuple


def merge_field_value(
    data_dict: Dict[str, Any],
    field: str,
    value: Any,
    method_name: str,
    providers: Dict[str, List[str]],
    context: Optional[Dict[str, Any]] = None,
) -> None:
    """Intelligently merge a field value into a data dictionary.

    Args:
        data_dict: Dictionary to update with the field value
        field: Field name to update
        value: New value to merge
        method_name: Name of the method providing the value
        providers: Dictionary tracking which methods provide which fields
        context: Optional context for error messages (e.g., gene name, site number)

    This function implements the following merge logic:
    1. If the field doesn't exist yet, simply add it
    2. If values are the same, no action needed
    3. If one value is 'NA', use the non-NA value
    4. If values differ and neither is 'NA', print a warning and use the latest value
    """
    # Track which methods provide this field
    if field not in providers:
        providers[field] = []
    if method_name not in providers[field]:
        providers[field].append(method_name)

    # If field doesn't exist yet, simply add it
    if field not in data_dict:
        data_dict[field] = value
        return

    # Field already exists, need to merge
    existing_value = data_dict[field]

    # Case 1: Values are the same - no action needed
    if existing_value == value:
        return

    # Case 2: One value is 'NA' - use the non-NA value
    if existing_value == "NA" and value != "NA":
        data_dict[field] = value
        return
    elif value == "NA" and existing_value != "NA":
        # Keep existing value
        return

    # Case 3: Both values differ and neither is 'NA' - print warning but keep latest
    context_str = ""
    if context:
        if "gene" in context:
            context_str += f" in gene {context['gene']}"
        if "site" in context:
            context_str += f" at site {context['site']}"

    print(f"WARNING: Conflicting values for field '{field}'{context_str}:")
    print(
        f"  Value 1: {existing_value} (from {', '.join(providers[field][:-1]) if len(providers[field]) > 1 else providers[field][0] if providers[field] else 'unknown'})"
    )
    print(f"  Value 2: {value} (from {method_name})")
    print(f"  Using value from {method_name}")

    # Update with the new value
    data_dict[field] = value


def merge_method_data(
    target_dict: Dict[str, Any],
    method_data: Dict[str, Any],
    method_name: str,
    providers: Dict[str, List[str]],
    context: Optional[Dict[str, Any]] = None,
) -> None:
    """Merge data from a method into a target dictionary.

    Args:
        target_dict: Dictionary to update with method data
        method_data: Data from the method to merge
        method_name: Name of the method providing the data
        providers: Dictionary tracking which methods provide which fields
        context: Optional context for error messages (e.g., gene name, site number)
    """
    for field, value in method_data.items():
        merge_field_value(target_dict, field, value, method_name, providers, context)


def ensure_ordered_fields(
    all_fields: Set[str], priority_fields: List[str]
) -> List[str]:
    """Ensure priority fields come first in the field list.

    Args:
        all_fields: Set of all fields to include
        priority_fields: List of fields that should come first (in order)

    Returns:
        Ordered list of fields with priority fields first
    """
    ordered_fields = []

    # Add priority fields first (in the specified order)
    for field in priority_fields:
        if field in all_fields:
            ordered_fields.append(field)

    # Add remaining fields
    for field in sorted(all_fields):
        if field not in priority_fields:
            ordered_fields.append(field)

    return ordered_fields


def collect_method_fields(
    methods: List[Any],
    method_results: Dict[str, Dict[str, Any]],
    field_getter_name: str,
    base_fields: Set[str] = None,
) -> Set[str]:
    """Collect fields from methods that have results.

    Args:
        methods: List of method objects
        method_results: Dictionary of method results keyed by method name
        field_getter_name: Name of the method function that returns fields (e.g., 'get_summary_fields')
        base_fields: Set of base fields to always include

    Returns:
        Set of collected fields from methods with results
    """
    # Initialize with base fields if provided
    collected_fields = set(base_fields) if base_fields else set()

    # Collect fields from methods that have results
    for method in methods:
        if method.name in method_results and hasattr(method, field_getter_name):
            field_getter = getattr(method, field_getter_name)
            # Pass the results if the method supports it
            if "results" in field_getter.__code__.co_varnames:
                collected_fields.update(field_getter(method_results[method.name]))
            else:
                collected_fields.update(field_getter())

    return collected_fields


def validate_fields(
    expected_fields: Set[str],
    output_fields: Set[str],
    context: str = None,
    entity_name: str = None,
) -> None:
    """Validate that all expected fields are present in the output fields.

    Args:
        expected_fields: Set of fields that are expected
        output_fields: Set of fields that are actually in the output
        context: Context for the validation (e.g., 'summary', 'site')
        entity_name: Name of the entity being validated (e.g., gene name)
    """
    if not expected_fields.issubset(output_fields):
        missing_fields = expected_fields - output_fields
        context_str = f" {context}" if context else ""
        entity_str = f" for {entity_name}" if entity_name else ""
        print(f"Warning: Missing{context_str} fields{entity_str}: {missing_fields}")


def detect_comparison_groups(
    method_results: Dict[str, Dict[str, Any]],
    default_groups: Optional[List[str]] = None,
    methods_to_check: Optional[List[str]] = None,
    gene_name: Optional[str] = None,
) -> Tuple[List[str], Dict[str, List[str]]]:
    """Detect comparison groups from method results.

    Args:
        method_results: Dictionary of method results keyed by method name
        default_groups: Default comparison groups to use if none are detected
        methods_to_check: List of method names to check for comparison groups.
                         If None, defaults to ['CFEL', 'RELAX']

    Returns:
        Tuple containing:
            - List of detected comparison groups
            - Dictionary with detected groups by method

    Raises:
        ValueError: If groups from different methods are inconsistent
    """
    # Default to CFEL and RELAX if no methods specified
    if methods_to_check is None:
        methods_to_check = ["CFEL", "RELAX"]

    # Initialize storage for detected groups
    detected_by_method: Dict[str, List[str]] = {}

    # Method-specific detection functions
    def detect_cfel_groups() -> List[str]:
        """Extract comparison groups from CFEL results."""
        if "CFEL" in method_results and method_results["CFEL"].get("tested", {}).get(
            "0"
        ):
            cfel_tested = method_results["CFEL"]["tested"]["0"]
            detected_groups: Set[str] = set()

            # Extract unique group tags from tested branches
            for branch, tag in cfel_tested.items():
                detected_groups.add(tag)

            if detected_groups:
                groups = list(detected_groups)
                gene_info = f" for {gene_name}" if gene_name else ""
                print(
                    f"Detected comparison groups from CFEL results{gene_info}: {', '.join(groups)}"
                )
                return groups
        return []

    def detect_relax_groups() -> List[str]:
        """Extract comparison groups from RELAX results."""
        if "RELAX" in method_results and method_results["RELAX"].get(
            "test results", {}
        ).get("relaxation or intensification parameter"):
            relax_k = method_results["RELAX"]["test results"][
                "relaxation or intensification parameter"
            ]
            if isinstance(relax_k, dict):
                detected_groups = [k for k in relax_k.keys() if k != "overall"]
                if detected_groups:
                    gene_info = f" for {gene_name}" if gene_name else ""
                    print(
                        f"Detected comparison groups from RELAX results{gene_info}: {', '.join(detected_groups)}"
                    )
                    return detected_groups
        return []

    # Map of method names to their detection functions
    detection_functions = {
        "CFEL": detect_cfel_groups,
        "RELAX": detect_relax_groups,
        # Add more methods here as they are implemented
    }

    # Detect groups from each specified method
    for method in methods_to_check:
        if method in detection_functions:
            groups = detection_functions[method]()
            if groups:
                detected_by_method[method] = groups
        else:
            print(f"Warning: No detection function available for method '{method}'")

    # Determine which groups to use
    all_detected_groups = list(detected_by_method.values())

    if not all_detected_groups:
        # No groups detected, use defaults
        if default_groups:
            comparison_groups = default_groups.copy()
            gene_info = f" for {gene_name}" if gene_name else ""
            print(
                f"Using default comparison groups{gene_info}: {', '.join(default_groups)}"
            )
        else:
            comparison_groups = []
        return comparison_groups, detected_by_method

    # If only one method detected groups, use those
    if len(all_detected_groups) == 1:
        return all_detected_groups[0], detected_by_method

    # Multiple methods detected groups, check consistency
    first_groups = all_detected_groups[0]
    first_method = list(detected_by_method.keys())[0]

    for i, groups in enumerate(all_detected_groups[1:], 1):
        method = list(detected_by_method.keys())[i]

        # Check if groups are consistent (ignoring order)
        if set(first_groups) != set(groups):
            # Groups are inconsistent, raise error
            error_msg = "ERROR: Inconsistent comparison groups between methods\n"
            for m, g in detected_by_method.items():
                error_msg += f"  {m} groups: {', '.join(g)}\n"
            print(error_msg)
            raise ValueError(error_msg)
        else:
            gene_info = f" for {gene_name}" if gene_name else ""
            print(
                f"Comparison groups are consistent{gene_info} between {first_method} and {method} methods"
            )

    # All groups are consistent, use the first one (arbitrary choice)
    return first_groups, detected_by_method
