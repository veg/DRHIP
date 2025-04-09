"""
Tests for HyPhy analysis methods.
"""

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

def test_busted_method_processing(mock_busted_results):
    """Test BUSTED method result processing."""
    method = BustedMethod()
    results = method.process_results(mock_busted_results)
    
    assert results['busted_pvalue'] == 0.01
    assert results['busted_lrt'] == 10.5
    assert results['busted_evidence']
    assert results['busted_omega3'] == 3.5
    assert results['busted_omega3_weight'] == 0.1
    assert results['total_branch_length'] == 0.3
    assert results['tree_length_ratio'] == 1.0

def test_relax_method_processing(mock_relax_results):
    """Test RELAX method result processing."""
    method = RelaxMethod()
    results = method.process_results(mock_relax_results)
    
    assert results['relax_pvalue'] == 0.05
    assert results['relax_k'] == 1.5
    assert results['relax_lrt'] == 8.2
    assert results['relax_k_significant']

def test_cfel_method_processing(mock_cfel_results):
    """Test CFEL method result processing."""
    method = CfelMethod()
    results = method.process_results(mock_cfel_results)
    
    assert 'by_type' in results
    assert 'data_rows' in results
    assert 'beta_idx_map' in results
    assert 'subs_idx_map' in results
    
    # Test site data processing
    site_data = method.process_site_data(results, ['clade1', 'clade2'])
    assert len(site_data) == 2  # Two sites in mock data
    
    site1_data = site_data[1]
    assert site1_data['beta_clade1'] == 0.5
    assert site1_data['subs_clade1'] == 2
    assert site1_data['beta_clade2'] == 1.0
    assert site1_data['subs_clade2'] == 0

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

def test_method_file_paths(temp_dir):
    """Test that methods generate correct file paths."""
    methods = [
        BustedMethod(),
        RelaxMethod(),
        CfelMethod(),
        FelMethod(),
        MemeMethod(),
        PrimeMethod()
    ]
    
    gene = 'test_gene'
    for method in methods:
        path = method.get_file_path(temp_dir, gene)
        assert gene in path
        assert method.file_suffix in path
