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
