"""
Tests for HyPhy analysis methods.
"""

import os

from hyphy_results_toolkit.methods import (
    HyPhyMethodRegistry,
    BustedMethod,
    RelaxMethod,
    CfelMethod,
    FelMethod,
    MemeMethod,
    PrimeMethod
)

def test_method_registry_initialization():
    """Test that the method registry initializes with all methods."""
    registry = HyPhyMethodRegistry()
    methods = registry.get_all_methods()
    
    # Check that all methods are registered
    method_names = {m.name for m in methods}
    expected_names = {'BUSTED', 'RELAX', 'CFEL', 'FEL', 'MEME', 'PRIME'}
    assert method_names == expected_names

def test_busted_method_processing(real_busted_results):
    """Test BUSTED method result processing."""
    method = BustedMethod()
    results = method.process_results(real_busted_results)
    
    # Basic validation of results
    assert 'BUSTED_pval' in results
    assert 'BUSTED_omega3' in results
    assert 'BUSTED_prop_sites_in_omega3' in results
    
    # Check data types
    assert isinstance(results['BUSTED_pval'], (float, int)) or results['BUSTED_pval'] == 'NA'

def test_fel_method_processing(real_fel_results):
    """Test FEL method result processing."""
    method = FelMethod()
    results = method.process_results(real_fel_results)
    
    # Basic validation of results - check for current field names
    assert 'N' in results
    assert 'positive_sites' in results
    assert 'negative_sites' in results
    
    # Test site data processing
    site_data = method.process_site_data(real_fel_results)
    
    # Check site data structure
    if len(site_data) > 0:
        first_site = min(site_data.keys())
        site_info = site_data[first_site]
        
        # Just verify that we have some data for the site
        assert isinstance(site_info, dict)
        assert len(site_info) > 0

def test_meme_method_processing(real_meme_results):
    """Test MEME method result processing."""
    method = MemeMethod()
    results = method.process_results(real_meme_results)
    
    # Verify the method returns a dictionary
    assert isinstance(results, dict)
    
    # Test site data processing
    site_data = method.process_site_data(real_meme_results)
    
    # Check site data structure if available
    if len(site_data) > 0:
        first_site = min(site_data.keys())
        site_info = site_data[first_site]
        
        # Just verify that we have some data for the site
        assert isinstance(site_info, dict)
        assert len(site_info) > 0

def test_method_field_generation():
    """Test that methods correctly generate field names."""
    registry = HyPhyMethodRegistry()
    
    # Test summary fields
    summary_fields = registry.get_all_summary_fields()
    
    # Check for presence of some fields from each method
    # Use the current field naming conventions
    assert 'BUSTED_omega3' in summary_fields or 'BUSTED_pval' in summary_fields
    assert 'positive_sites' in summary_fields or 'negative_sites' in summary_fields
    
    # Test site fields without comparison groups
    site_fields = registry.get_all_site_fields()
    
    # Check for some expected fields
    assert isinstance(site_fields, list)
    assert len(site_fields) > 0

def test_prime_method_processing(real_prime_results):
    """Test PRIME method result processing."""
    method = PrimeMethod()
    results = method.process_results(real_prime_results)
    
    # Verify the method returns a dictionary
    assert isinstance(results, dict)
    
    # Test site data processing
    site_data = method.process_site_data(real_prime_results)
    
    # Verify site data is a dictionary
    assert isinstance(site_data, dict)

def test_method_file_paths(results_dir):
    """Test that methods generate correct file paths."""
    methods = [
        BustedMethod(),
        RelaxMethod(),
        CfelMethod(),
        FelMethod(),
        MemeMethod(),
        PrimeMethod()
    ]
    
    gene = 'capsid_protein_C'
    for method in methods:
        if method.name in ['BUSTED', 'FEL', 'MEME', 'PRIME']:
            path = method.get_file_path(results_dir, gene)
            expected_file = f"{gene}.{method.file_suffix}"
            assert os.path.exists(path)
            assert expected_file in path

'''
# Mock data for RELAX tests
@pytest.fixture
def mock_relax_results_with_data():
    """Create mock RELAX results with valid data."""
    return {
        "test results": {
            "p-value": 0.01,
            "relaxation or intensification parameter": 0.5,
            "evidence": "relaxation",
            "test statistic": 8.5
        },
        "fits": {
            "Global MG94xREV": {
                "Rate Distributions": {
                    "foreground": [[0.5, 1.0]],
                    "background": [[0.8, 1.0]]
                }
            }
        },
        "branch attributes": {
            "0": {
                "branch1": {"Global MG94xREV": 0.1},
                "branch2": {"Global MG94xREV": 0.2},
                "branch3": {"Global MG94xREV": 0.3}
            }
        },
        "tested": {
            "0": {
                "branch1": "foreground",
                "branch2": "foreground",
                "branch3": "background"
            }
        }
    }

@pytest.fixture
def mock_relax_results_empty():
    """Create mock RELAX results with empty data."""
    return {}

def test_relax_method_with_comparison_groups(mock_relax_results_with_data):
    """Test RELAX method with comparison groups."""
    method = RelaxMethod()
    
    # Set comparison groups
    comparison_groups = ["foreground", "background"]
    method._comparison_groups = comparison_groups
    
    # Process comparison data
    results = method.process_comparison_data(mock_relax_results_with_data)
    
    # Verify results
    assert 'RELAX_overall_pval' in results
    assert 'RELAX_K' in results
    assert results['RELAX_overall_pval'] == 0.01
    assert results['RELAX_K'] == 0.5
    
    # Check comparison group fields
    assert method.get_comparison_group_summary_fields() == [
        'RELAX_overall_pval',
        'RELAX_K'
    ]
    
    # Verify no site data is produced
    site_data = method.process_site_data(mock_relax_results_with_data)
    assert site_data == {}

def test_relax_method_without_comparison_groups(mock_relax_results_with_data):
    """Test RELAX method without comparison groups."""
    method = RelaxMethod()
    
    # Process results without setting comparison groups
    with pytest.raises(ValueError, match="Comparison groups must be set"):
        method.process_comparison_data(mock_relax_results_with_data)
    
    # Verify standard results processing
    results = method.process_results(mock_relax_results_with_data)
    assert results == {}  # RELAX only uses comparison group fields

def test_relax_method_with_empty_data(mock_relax_results_empty):
    """Test RELAX method with empty data."""
    method = RelaxMethod()
    
    # Set comparison groups
    comparison_groups = ["foreground", "background"]
    method._comparison_groups = comparison_groups
    
    # Process comparison data with empty results
    results = method.process_comparison_data(mock_relax_results_empty)
    
    # Verify NA values for missing data
    assert results['RELAX_overall_pval'] == 'NA'
    assert results['RELAX_K'] == 'NA'

# Mock data for CFEL tests
@pytest.fixture
def mock_cfel_results_with_data():
    """Create mock CFEL results with valid data."""
    return {
        "tested": {
            "0": {
                "branch1": "foreground",
                "branch2": "foreground",
                "branch3": "background"
            }
        },
        "branch attributes": {
            "0": {
                "branch1": {"Global MG94xREV": 0.1},
                "branch2": {"Global MG94xREV": 0.2},
                "branch3": {"Global MG94xREV": 0.3}
            }
        },
        "fits": {
            "Global MG94xREV": {
                "Rate Distributions": {
                    "foreground": [[0.5, 1.0]],
                    "background": [[0.8, 1.0]]
                }
            }
        },
        "MLE": {
            "headers": [
                "Site",
                "α",
                "β (foreground)",
                "β (background)",
                "substitutions (foreground)",
                "substitutions (background)"
            ],
            "content": {
                "0": [
                    {"value": ["1", "0.1", "0.0", "0.5", "1", "0"]},
                    {"value": ["2", "0.2", "0.5", "0.0", "0", "2"]},
                    {"value": ["3", "0.3", "0.0", "0.0", "1", "1"]}
                ]
            }
        }
    }

@pytest.fixture
def mock_cfel_results_empty():
    """Create mock CFEL results with empty data."""
    return {}

def test_cfel_method_with_comparison_groups(mock_cfel_results_with_data):
    """Test CFEL method with comparison groups."""
    method = CfelMethod()
    
    # Set comparison groups
    comparison_groups = ["foreground", "background"]
    method._comparison_groups = comparison_groups
    
    # Process comparison data
    results = method.process_comparison_data(mock_cfel_results_with_data)
    
    # Verify results
    assert "foreground" in results
    assert "background" in results
    assert "group_N" in results["foreground"]
    assert "group_T" in results["foreground"]
    assert "group_dN/dS" in results["foreground"]
    
    # Check comparison group fields
    assert set(method.get_comparison_group_summary_fields()) == {
        'group_N',
        'group_T',
        'group_dN/dS',
        'group_nt_conserved',
        'group_aa_conserved'
    }
    
    # Test site-specific data
    site_data = method.process_comparison_site_data(mock_cfel_results_with_data)
    assert len(site_data) > 0
    
    # Check site fields
    assert set(method.get_comparison_group_site_fields()) == {'cfel_marker'}

def test_cfel_method_without_comparison_groups(mock_cfel_results_with_data):
    """Test CFEL method without comparison groups."""
    method = CfelMethod()
    
    # Process results without setting comparison groups
    with pytest.raises(ValueError, match="Comparison groups must be set"):
        method.process_comparison_data(mock_cfel_results_with_data)
    
    # Verify standard results processing
    results = method.process_results(mock_cfel_results_with_data)
    assert results == {}  # CFEL only uses comparison group fields
    
    # Verify site data processing without comparison groups
    site_data = method.process_site_data(mock_cfel_results_with_data)
    assert site_data == {}  # CFEL only uses comparison group site fields

def test_cfel_method_with_empty_data(mock_cfel_results_empty):
    """Test CFEL method with empty data."""
    method = CfelMethod()
    
    # Set comparison groups
    comparison_groups = ["foreground", "background"]
    method._comparison_groups = comparison_groups
    
    # Process comparison data with empty results
    results = method.process_comparison_data(mock_cfel_results_empty)
    
    # Verify empty results for missing data
    assert results == {}
'''
