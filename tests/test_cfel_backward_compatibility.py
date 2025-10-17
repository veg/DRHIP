"""
Tests for CFEL backward compatibility with older HyPhy versions.

These tests ensure that CFEL processing works correctly with older HyPhy outputs
that do not contain the 'substitutions' field, which was added in newer versions.
This prevents accidental breakage when processing legacy data.
"""

import json
import os
import warnings
from typing import Any, Dict

from drhip.methods.cfel import CfelMethod


def _load_cfel_results(
    comparison_results_dir: str, gene_basename: str
) -> Dict[str, Any]:
    """Load CFEL results from test data directory."""
    # Handle both with and without .CONTRASTFEL suffix in basename
    if gene_basename.endswith(".CONTRASTFEL") or gene_basename.endswith(
        ".CONTRASTFEL.no-subs"
    ):
        filename = f"{gene_basename}.json"
    else:
        filename = f"{gene_basename}.CONTRASTFEL.json"

    path = os.path.join(
        comparison_results_dir,
        "CONTRASTFEL",
        filename,
    )
    with open(path) as f:
        return json.load(f)


def test_cfel_validation_passes_without_substitutions(comparison_results_dir: str):
    """
    Test that CFEL validation passes for files without substitutions field.

    This ensures backward compatibility - older HyPhy versions don't have
    substitutions, and we should still be able to process those files.
    """
    # Load a file without substitutions
    gene_basename = "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL.no-subs"
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    # Verify it doesn't have substitutions
    assert (
        "substitutions" not in results
    ), "Test data should not have substitutions field"

    # Validation should pass
    method = CfelMethod()
    missing = method.validate_input_json(results)

    # Should have no missing required fields
    assert len(missing) == 0, f"Validation failed with missing fields: {missing}"


def test_cfel_process_results_without_substitutions(comparison_results_dir: str):
    """
    Test that process_results works correctly without substitutions field.

    The summary-level processing (diff_sites count) should work regardless
    of whether substitutions are present.
    """
    gene_basename = "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL.no-subs"
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    assert "substitutions" not in results

    method = CfelMethod()
    processed = method.process_results(results)

    # Should still calculate diff_sites
    assert "diff_sites" in processed
    assert isinstance(processed["diff_sites"], int)
    assert processed["diff_sites"] >= 0


def test_cfel_process_comparison_data_without_substitutions(
    comparison_results_dir: str,
):
    """
    Test that process_comparison_data works without substitutions field.

    Group-level summary data (N, T, dN/dS, aa_conserved) should be calculated
    independently of substitutions data.
    """
    gene_basename = "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL.no-subs"
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    assert "substitutions" not in results

    # Get comparison groups
    tested = results.get("tested", {}).get("0", {})
    groups = sorted(set(tested.values()))

    method = CfelMethod().set_comparison_groups(groups)
    comp_data = method.process_comparison_data(results)

    # Should have data for each group
    for group in groups:
        assert group in comp_data
        gdata = comp_data[group]

        # All required fields should be present
        assert "group_N" in gdata
        assert "group_T" in gdata
        assert "group_dN/dS" in gdata
        assert "group_aa_conserved" in gdata

        # Values should be valid
        assert isinstance(gdata["group_N"], int)
        assert gdata["group_N"] > 0
        assert isinstance(gdata["group_T"], (int, float))
        assert gdata["group_T"] >= 0


def test_cfel_process_comparison_site_data_without_substitutions_warns(
    comparison_results_dir: str,
):
    """
    Test that process_comparison_site_data emits a warning when substitutions are missing.

    This is critical for backward compatibility - users should be warned that
    composition/substitutions/majority_residue fields won't be populated.
    """
    gene_basename = "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL.no-subs"
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    assert "substitutions" not in results

    tested = results.get("tested", {}).get("0", {})
    groups = sorted(set(tested.values()))

    method = CfelMethod().set_comparison_groups(groups)

    # Should emit a UserWarning about missing substitutions
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _ = method.process_comparison_site_data(results)

        # Check that a warning was issued
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert "substitutions" in str(w[0].message).lower()
        assert "older HyPhy versions" in str(w[0].message)


def test_cfel_process_comparison_site_data_without_substitutions_has_core_fields(
    comparison_results_dir: str,
):
    """
    Test that core site fields (marker, beta) are still populated without substitutions.

    Even without substitutions, we should still get cfel_marker and cfel_beta values
    since those come from the MLE content, not the substitutions field.
    """
    gene_basename = "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL.no-subs"
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    assert "substitutions" not in results

    tested = results.get("tested", {}).get("0", {})
    groups = sorted(set(tested.values()))

    method = CfelMethod().set_comparison_groups(groups)

    # Suppress the expected warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        site_data = method.process_comparison_site_data(results)

    # Should still have site data
    assert isinstance(site_data, dict)
    assert len(site_data) > 0

    # Check first few sites
    for site_id in list(site_data.keys())[:5]:
        for group in groups:
            assert group in site_data[site_id]
            gdata = site_data[site_id][group]

            # Core fields from MLE should be present
            assert "cfel_marker" in gdata
            assert "cfel_beta" in gdata

            # These fields should NOT be present (they come from substitutions)
            assert "composition" not in gdata
            assert "substitutions" not in gdata
            assert "majority_residue" not in gdata


def test_cfel_comparison_with_and_without_substitutions(comparison_results_dir: str):
    """
    Test that files with and without substitutions can both be processed.

    This ensures that the presence/absence of substitutions doesn't break processing.
    Note: The test files may be from different HyPhy versions with different column
    layouts, so we don't compare exact values, just that both process successfully.
    """
    # Load both versions of the same gene
    gene_base = (
        "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL"
    )

    results_with_subs = _load_cfel_results(comparison_results_dir, gene_base)
    results_without_subs = _load_cfel_results(
        comparison_results_dir, f"{gene_base}.no-subs"
    )

    assert "substitutions" in results_with_subs
    assert "substitutions" not in results_without_subs

    # Get comparison groups (should be the same)
    tested_with = results_with_subs.get("tested", {}).get("0", {})
    groups = sorted(set(tested_with.values()))

    # Process both
    method_with = CfelMethod().set_comparison_groups(groups)
    method_without = CfelMethod().set_comparison_groups(groups)

    # Both should process successfully
    summary_with = method_with.process_results(results_with_subs)
    summary_without = method_without.process_results(results_without_subs)

    # Both should have diff_sites
    assert "diff_sites" in summary_with
    assert "diff_sites" in summary_without
    assert isinstance(summary_with["diff_sites"], int)
    assert isinstance(summary_without["diff_sites"], int)

    # Compare site data (suppress warning for no-subs version)
    site_data_with = method_with.process_comparison_site_data(results_with_subs)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        site_data_without = method_without.process_comparison_site_data(
            results_without_subs
        )

    # Both should have site data
    assert len(site_data_with) > 0
    assert len(site_data_without) > 0

    # Check that 'with' version has substitution-derived fields
    first_site_with = list(site_data_with.keys())[0]
    first_group = groups[0]
    data_with = site_data_with[first_site_with][first_group]

    assert "cfel_marker" in data_with
    assert "cfel_beta" in data_with
    assert "composition" in data_with
    assert "substitutions" in data_with
    assert "majority_residue" in data_with

    # Check that 'without' version lacks substitution-derived fields
    first_site_without = list(site_data_without.keys())[0]
    data_without = site_data_without[first_site_without][first_group]

    assert "cfel_marker" in data_without
    assert "cfel_beta" in data_without
    assert "composition" not in data_without
    assert "substitutions" not in data_without
    assert "majority_residue" not in data_without


def test_cfel_no_regression_on_newer_files(comparison_results_dir: str):
    """
    Test that files WITH substitutions still work correctly.

    This is a regression test to ensure our backward compatibility changes
    don't break processing of newer HyPhy outputs.
    """
    gene_basename = (
        "pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1.CONTRASTFEL"
    )
    results = _load_cfel_results(comparison_results_dir, gene_basename)

    # Verify it HAS substitutions
    assert "substitutions" in results, "Test data should have substitutions field"

    tested = results.get("tested", {}).get("0", {})
    groups = sorted(set(tested.values()))

    method = CfelMethod().set_comparison_groups(groups)

    # Should NOT emit a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        site_data = method.process_comparison_site_data(results)

        # No warnings should be issued
        assert len(w) == 0, f"Unexpected warnings: {[str(x.message) for x in w]}"

    # Should have all fields including substitution-derived ones
    assert len(site_data) > 0

    for site_id in list(site_data.keys())[:5]:
        for group in groups:
            gdata = site_data[site_id][group]

            # All fields should be present
            assert "cfel_marker" in gdata
            assert "cfel_beta" in gdata
            assert "composition" in gdata
            assert "substitutions" in gdata
            assert "majority_residue" in gdata
