"""
Tests covering PR #11 for CFEL comparison outputs:
- process_comparison_data: group_N, group_T, group_dN/dS, group_aa_conserved
- process_comparison_site_data: composition, substitutions, majority_residue presence
"""

import json
import os
from typing import Dict, Any, List

import math
import pytest

from drhip.methods.cfel import CfelMethod


def _load_cfel_results(comparison_results_dir: str, gene_basename: str) -> Dict[str, Any]:
    path = os.path.join(
        comparison_results_dir,
        "CONTRASTFEL",
        f"{gene_basename}.CONTRASTFEL.json",
    )
    with open(path, "r") as f:
        return json.load(f)


def _list_cfel_basenames(comparison_results_dir: str) -> List[str]:
    dir_path = os.path.join(comparison_results_dir, "CONTRASTFEL")
    basenames = []
    if os.path.isdir(dir_path):
        for file in os.listdir(dir_path):
            if file.endswith(".CONTRASTFEL.json"):
                basenames.append(file[: -len(".CONTRASTFEL.json")])
    return sorted(basenames)


def _get_header_indices(results: Dict[str, Any]) -> Dict[str, int]:
    header_indices: Dict[str, int] = {}
    for i, (name, _) in enumerate(results["MLE"]["headers"]):
        header_indices[name] = i
    return header_indices


def _expected_group_counts_and_T(results: Dict[str, Any], groups: List[str]):
    tested = results.get("tested", {}).get("0", {})
    by_type = {g: [] for g in groups}
    for branch, tag in tested.items():
        if tag in by_type:
            by_type[tag].append(branch)

    # Compute expected N and T
    branch_attrs = results.get("branch attributes", {}).get("0", {})
    expected = {}
    for g, branches in by_type.items():
        T = 0.0
        for bn in branches:
            try:
                T += float(branch_attrs[bn]["Global MG94xREV"])  # type: ignore[index]
            except Exception:
                T += 0.0
        expected[g] = {"N": len(branches), "T": T}
    return expected


def _expected_group_dNdS(results: Dict[str, Any], groups: List[str]):
    rd = results.get("fits", {}).get("Global MG94xREV", {}).get("Rate Distributions", {})
    out = {}
    for g in groups:
        found = None
        for k, r in rd.items():
            if g in k and isinstance(r, list) and r and isinstance(r[0], list) and r[0]:
                found = float(r[0][0])
                break
        out[g] = found if found is not None else "NA"
    return out


def _expected_group_aa_conserved(results: Dict[str, Any], groups: List[str]):
    headers = _get_header_indices(results)
    rows = results.get("MLE", {}).get("content", {}).get("0", [])
    # Map group -> beta column index
    beta_idx_map: Dict[str, int] = {}
    for k, idx in headers.items():
        lk = k.lower()
        if lk.startswith("beta ") and "(" in k and ")" in k:
            group_name = k[k.find("(") + 1 : k.find(")")]
            if group_name in groups:
                beta_idx_map[group_name] = idx

    out = {}
    for g in groups:
        idx = beta_idx_map.get(g, -1)
        if idx < 0:
            out[g] = "NA"
            continue
        count = 0
        for row in rows:
            if idx < len(row):
                try:
                    beta_val = float(row[idx])
                    if beta_val == 0.0:
                        count += 1
                except (ValueError, TypeError):
                    continue
        out[g] = count
    return out


def test_cfel_process_comparison_data_fields_and_values(comparison_results_dir: str):
    basenames = _list_cfel_basenames(comparison_results_dir)
    assert basenames, "No CFEL comparison datasets found for testing"

    for gene_basename in basenames:
        results = _load_cfel_results(comparison_results_dir, gene_basename)

        # Determine comparison groups from tested map
        tested = results.get("tested", {}).get("0", {})
        groups = sorted(set(tested.values()))
        assert groups, "No comparison groups detected in CFEL test data"

        method = CfelMethod().set_comparison_groups(groups)
        comp_data = method.process_comparison_data(results)

        # Must have an entry per group
        for g in groups:
            assert g in comp_data, f"Missing comparison group {g} in summary data"
            gdata = comp_data[g]
            # Required keys
            for key in ["group_N", "group_T", "group_dN/dS", "group_aa_conserved"]:
                assert key in gdata, f"Missing {key} for group {g}"

        # Validate N and T
        expected_counts = _expected_group_counts_and_T(results, groups)
        for g in groups:
            assert comp_data[g]["group_N"] == expected_counts[g]["N"]
            # Floating point tolerance for T
            assert math.isclose(float(comp_data[g]["group_T"]), expected_counts[g]["T"], rel_tol=1e-9, abs_tol=1e-9)

        # Validate dN/dS
        expected_dNdS = _expected_group_dNdS(results, groups)
        for g in groups:
            if expected_dNdS[g] == "NA":
                assert comp_data[g]["group_dN/dS"] == "NA"
            else:
                assert math.isclose(float(comp_data[g]["group_dN/dS"]), float(expected_dNdS[g]), rel_tol=1e-9, abs_tol=1e-9)

        # Validate AA conserved count
        expected_aa = _expected_group_aa_conserved(results, groups)
        for g in groups:
            if expected_aa[g] == "NA":
                assert comp_data[g]["group_aa_conserved"] == "NA"
            else:
                assert comp_data[g]["group_aa_conserved"] == expected_aa[g], f"Comparison group {g} doesn't match expected AA conserved count"


def test_cfel_process_comparison_site_data_has_composition_and_substitutions(comparison_results_dir: str):
    basenames = _list_cfel_basenames(comparison_results_dir)
    assert basenames, "No CFEL comparison datasets found for testing"

    # Only sample first dataset to keep runtime sane
    gene_basename = basenames[0]
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    tested = results.get("tested", {}).get("0", {})
    groups = sorted(set(tested.values()))
    method = CfelMethod().set_comparison_groups(groups)

    site_data = method.process_comparison_site_data(results)
    assert isinstance(site_data, dict)
    assert site_data, "No site data returned"

    # Inspect first 30 rows (first 10 sites, with 3 groups per site) present in MLE content or substitutions
    rows = results.get("MLE", {}).get("content", {}).get("0", [])
    num_sites = min(30, len(rows))
    for i in range(num_sites):
        site_id = str(i + 1)
        assert site_id in site_data
        for g in groups:
            assert g in site_data[site_id]
            gdata = site_data[site_id][g]
            # Presence of keys
            for key in ["cfel_marker", "cfel_beta", "composition", "substitutions", "majority_residue"]:
                # currently throwing errors when a site has no subs & therefore comp/subs/major residue isn't calculated 
                assert key in gdata, f"Missing key {key} for group {g} at site {site_id}"
            # Basic format checks
            assert isinstance(gdata["composition"], str)
            assert isinstance(gdata["substitutions"], str)
            assert isinstance(gdata["majority_residue"], str)
            assert gdata["majority_residue"] == '-' or len(gdata["majority_residue"]) == 1
