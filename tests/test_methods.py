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
    assert 'busted_pvalue' in results
    assert 'busted_lrt' in results
    assert 'busted_evidence' in results
    assert 'busted_omega3' in results
    assert 'busted_omega3_weight' in results
    assert 'total_branch_length' in results
    assert 'tree_length_ratio' in results
    
    # Check evidence type
    assert isinstance(results['busted_evidence'], bool)

def test_fel_method_processing(real_fel_results):
    """Test FEL method result processing."""
    method = FelMethod()
    results = method.process_results(real_fel_results)
    
    # Basic validation of results
    assert 'fel_sites_tested' in results
    assert 'fel_sites_positive_selection' in results
    assert 'fel_sites_negative_selection' in results
    assert 'fel_version' in results
    assert 'fel_timestamp' in results
    
    # Test site data processing
    site_data = method.process_site_data(real_fel_results)
    assert len(site_data) > 0
    
    # Check a specific site
    first_site = min(site_data.keys())
    site_info = site_data[first_site]
    assert 'fel_alpha' in site_info
    assert 'fel_beta' in site_info
    assert 'fel_pvalue' in site_info
    assert 'fel_selection' in site_info
    assert site_info['fel_selection'] in ['positive', 'negative', 'neutral']

def test_meme_method_processing(real_meme_results):
    """Test MEME method result processing."""
    method = MemeMethod()
    results = method.process_results(real_meme_results)
    
    # Basic validation of results
    assert 'meme_sites_tested' in results
    assert 'meme_sites_selection' in results
    assert 'meme_version' in results
    assert 'meme_timestamp' in results
    
    # Test site data processing
    site_data = method.process_site_data(real_meme_results)
    assert len(site_data) > 0
    
    # Check a specific site
    first_site = min(site_data.keys())
    site_info = site_data[first_site]
    assert 'meme_alpha' in site_info
    assert 'meme_beta_neg' in site_info
    assert 'meme_beta_plus' in site_info
    assert 'meme_weight' in site_info
    assert 'meme_pvalue' in site_info
    assert 'meme_selection' in site_info

def test_method_field_generation():
    """Test that methods correctly generate field names."""
    registry = HyPhyMethodRegistry()
    
    # Test summary fields
    summary_fields = registry.get_all_summary_fields()
    assert 'busted_pvalue' in summary_fields
    assert 'relax_k' in summary_fields
    assert 'fel_sites_tested' in summary_fields
    assert 'meme_sites_selection' in summary_fields
    
    # Test site fields
    site_fields = registry.get_all_site_fields(['clade1', 'clade2'])
    assert 'beta_clade1' in site_fields
    assert 'subs_clade2' in site_fields
    assert 'fel_alpha' in site_fields
    assert 'meme_pvalue' in site_fields

def test_prime_method_processing(real_prime_results):
    """Test PRIME method result processing."""
    method = PrimeMethod()
    results = method.process_results(real_prime_results)
    
    # Basic validation of results
    assert 'prime_sites_tested' in results
    assert 'prime_version' in results
    assert 'prime_timestamp' in results
    
    # Check property-specific fields
    for prop in PrimeMethod.PROPERTIES:
        prop_key = prop.lower().replace(' ', '_')
        assert f'prime_{prop_key}_conserved' in results
        assert f'prime_{prop_key}_altered' in results
    
    # Test site data processing
    site_data = method.process_site_data(real_prime_results)
    assert len(site_data) > 0

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
