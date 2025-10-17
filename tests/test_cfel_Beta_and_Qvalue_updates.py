"""
Tests covering PR #7 (CFEL fixes):
- diff_sites counting from Q-value column
- cfel_marker formatting per site/group based on Q-value
- cfel_beta per-group value formatting
"""

import json
import os
import re
from typing import Any, Dict

import pytest

from drhip.methods.cfel import CfelMethod


def _load_cfel_results(
    comparison_results_dir: str, gene_basename: str
) -> Dict[str, Any]:
    path = os.path.join(
        comparison_results_dir,
        "CONTRASTFEL",
        f"{gene_basename}.CONTRASTFEL.json",
    )
    with open(path) as f:
        return json.load(f)


def _get_header_indices(results: Dict[str, Any]) -> Dict[str, int]:
    header_indices: Dict[str, int] = {}
    for i, (name, _) in enumerate(results["MLE"]["headers"]):
        header_indices[name] = i
    return header_indices


def _list_cfel_basenames(comparison_results_dir: str):
    dir_path = os.path.join(comparison_results_dir, "CONTRASTFEL")
    basenames = []
    if os.path.isdir(dir_path):
        for file in os.listdir(dir_path):
            if file.endswith(".CONTRASTFEL.json"):
                basenames.append(file[: -len(".CONTRASTFEL.json")])
    return sorted(basenames)


def test_cfel_diff_sites_count_from_qvalue(comparison_results_dir: str):
    """CfelMethod.process_results should count diff_sites as rows with Q-value <= 0.20."""
    basenames = _list_cfel_basenames(comparison_results_dir)
    assert basenames, "No CFEL comparison datasets found for testing"

    for gene_basename in basenames:
        results = _load_cfel_results(comparison_results_dir, gene_basename)

        # Compute expected diff_sites directly from raw JSON
        headers = _get_header_indices(results)
        q_idx = None
        # Q-value header may vary in case
        for k, idx in headers.items():
            if k.lower().startswith("q-value"):
                q_idx = idx
                break
        assert q_idx is not None, "Q-value column not found in CFEL headers"

        expected = 0
        for row in results["MLE"]["content"]["0"]:
            if q_idx < len(row):
                try:
                    q = float(row[q_idx])
                except (ValueError, TypeError):
                    continue
                if q <= 0.20:
                    expected += 1

        method = CfelMethod()
        processed = method.process_results(results)
        assert "diff_sites" in processed
        assert processed["diff_sites"] == expected


def test_cfel_marker_and_beta_formatting(comparison_results_dir: str):
    """
    CfelMethod.process_comparison_site_data should:
    - set cfel_marker to formatted Q-value ("{q:.3f}") when q <= 0.20, else '-'
    - set cfel_beta to per-group formatted value:
        * '0.000' for exact 0
        * scientific notation with 3 decimals if |beta|<1e-3 or >= 1e3
        * fixed-point with 4 decimals otherwise
    """
    basenames = _list_cfel_basenames(comparison_results_dir)
    assert basenames, "No CFEL comparison datasets found for testing"

    for gene_basename in basenames:
        results = _load_cfel_results(comparison_results_dir, gene_basename)

        # Determine comparison groups from tested map
        tested = results.get("tested", {}).get("0", {})
        groups = sorted(set(tested.values()))
        assert groups, "No comparison groups detected in CFEL test data"

        method = CfelMethod().set_comparison_groups(groups)
        site_data = method.process_comparison_site_data(results)

        # Build our own header index map for beta and Q-value
        headers = _get_header_indices(results)

        # Map group -> beta column index using headers like 'beta (group)'
        beta_idx_map: Dict[str, int] = {}
        for k, idx in headers.items():
            lk = k.lower()
            if lk.startswith("beta ") and "(" in k and ")" in k:
                group_name = k[k.find("(") + 1 : k.find(")")]
                if group_name in groups:
                    beta_idx_map[group_name] = idx

        # Locate the Q-value column index
        q_idx = None
        for k, idx in headers.items():
            if k.lower().startswith("q-value"):
                q_idx = idx
                break
        assert q_idx is not None, "Q-value column not found in CFEL headers"

        rows = results["MLE"]["content"]["0"]

        # Validate a sample of sites (e.g., first 30 or all if fewer)
        sample_count = min(30, len(rows))
        for site_zero_idx in range(sample_count):
            site_id = str(site_zero_idx + 1)
            row = rows[site_zero_idx]
            # Expected marker from Q-value
            expected_marker = "NA"
            if q_idx < len(row):
                try:
                    q = float(row[q_idx])
                    expected_marker = f"{q:.3f}" if q <= 0.20 else "-"
                except (ValueError, TypeError):
                    expected_marker = "NA"

            for group in groups:
                assert (
                    site_id in site_data
                ), f"Missing site {site_id} in comparison data"
                assert (
                    group in site_data[site_id]
                ), f"Missing group {group} at site {site_id}"

                gdata = site_data[site_id][group]

                # cfel_marker formatting
                assert "cfel_marker" in gdata
                marker = gdata["cfel_marker"]
                # Either '-', or a numeric string with exactly 3 decimals
                assert marker == expected_marker or marker == "NA"
                if marker not in ("-", "NA"):
                    assert re.fullmatch(r"\d*\.\d{3}", marker) is not None

                # cfel_beta formatting
                assert "cfel_beta" in gdata
                beta_str = gdata["cfel_beta"]
                # Compute expected formatting from raw value if available
                idx = beta_idx_map.get(group, -1)
                if 0 <= idx < len(row):
                    try:
                        beta_val = float(row[idx])
                        if beta_val == 0.0:
                            expected_beta = "0.000"
                        else:
                            ab = abs(beta_val)
                            if ab < 1e-3 or ab >= 1e3:
                                expected_beta = f"{beta_val:.3e}"
                            else:
                                expected_beta = f"{beta_val:.4f}"
                        assert beta_str == expected_beta
                    except (ValueError, TypeError):
                        assert beta_str == "NA"
                else:
                    assert beta_str == "NA"


def test_cfel_get_comparison_group_site_fields_contains_expected():
    fields = set(CfelMethod.get_comparison_group_site_fields())
    for f in [
        "cfel_marker",
        "cfel_beta",
        "composition",
        "substitutions",
        "majority_residue",
    ]:
        assert f in fields
