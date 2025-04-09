"""
Tests for process_gene module.
"""

import os
import csv
from hyphy_results_toolkit.parsers.process_gene import process_gene

def test_process_gene_output_files(temp_dir, mock_results_dir, write_mock_results):
    """Test that process_gene creates expected output files."""
    # Write mock results for test gene
    write_mock_results('test_gene')
    
    # Create output directory
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process gene
    process_gene('test_gene', mock_results_dir, output_dir, temp_dir)
    
    # Check that output files exist
    assert os.path.exists(os.path.join(output_dir, 'test_gene_summary.csv'))
    assert os.path.exists(os.path.join(output_dir, 'test_gene_sites.csv'))

def test_process_gene_summary_content(temp_dir, mock_results_dir, write_mock_results):
    """Test that process_gene produces correct summary content."""
    # Write mock results for test gene
    write_mock_results('test_gene')
    
    # Create output directory
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process gene
    process_gene('test_gene', mock_results_dir, output_dir, temp_dir)
    
    # Read summary file
    summary_file = os.path.join(output_dir, 'test_gene_summary.csv')
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f)
        row = next(reader)
        
        # Check BUSTED results
        assert float(row['busted_pvalue']) == 0.01
        assert float(row['busted_lrt']) == 10.5
        assert float(row['busted_omega3']) == 3.5
        
        # Check RELAX results
        assert float(row['relax_pvalue']) == 0.05
        assert float(row['relax_k']) == 1.5
        assert float(row['relax_lrt']) == 8.2

def test_process_gene_sites_content(temp_dir, mock_results_dir, write_mock_results):
    """Test that process_gene produces correct site-specific content."""
    # Write mock results for test gene
    write_mock_results('test_gene')
    
    # Create output directory
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process gene
    process_gene('test_gene', mock_results_dir, output_dir, temp_dir)
    
    # Read sites file
    sites_file = os.path.join(output_dir, 'test_gene_sites.csv')
    with open(sites_file, 'r') as f:
        reader = csv.DictReader(f)
        sites = list(reader)
        
        # Check first site data
        site1 = next(s for s in sites if s['site'] == '1')
        assert float(site1['beta_clade1']) == 0.5
        assert float(site1['subs_clade1']) == 2.0
        assert float(site1['beta_clade2']) == 1.0
        assert float(site1['subs_clade2']) == 0.0

def test_process_gene_missing_required_method(temp_dir, mock_results_dir, write_mock_results):
    """Test that process_gene handles missing required method results."""
    # Write incomplete results (missing BUSTED)
    gene = 'test_gene'
    
    # Write only RELAX and CFEL results
    with open(os.path.join(mock_results_dir, 'concat', 'RELAX', f'{gene}.RELAX.json'), 'w') as f:
        f.write('{}')
    with open(os.path.join(mock_results_dir, 'concat', 'contrastFEL', f'{gene}.CFEL.json'), 'w') as f:
        f.write('{}')
    
    # Create output directory
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process gene
    process_gene(gene, mock_results_dir, output_dir, temp_dir)
    
    # Check that output files were not created
    assert not os.path.exists(os.path.join(output_dir, f'{gene}_summary.csv'))
    assert not os.path.exists(os.path.join(output_dir, f'{gene}_sites.csv'))

def test_process_gene_thread_safety(temp_dir, mock_results_dir, write_mock_results):
    """Test that process_gene is thread-safe when writing output."""
    import threading
    
    # Write mock results for multiple genes
    genes = ['gene1', 'gene2', 'gene3']
    for gene in genes:
        write_mock_results(gene)
    
    # Create output directory
    output_dir = os.path.join(temp_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process genes in parallel
    threads = []
    for gene in genes:
        thread = threading.Thread(
            target=process_gene,
            args=(gene, mock_results_dir, output_dir, temp_dir)
        )
        threads.append(thread)
        thread.start()
    
    # Wait for all threads to complete
    for thread in threads:
        thread.join()
    
    # Check that all output files exist and are valid
    for gene in genes:
        summary_file = os.path.join(output_dir, f'{gene}_summary.csv')
        sites_file = os.path.join(output_dir, f'{gene}_sites.csv')
        
        assert os.path.exists(summary_file)
        assert os.path.exists(sites_file)
        
        # Verify files are valid CSV
        with open(summary_file, 'r') as f:
            reader = csv.DictReader(f)
            assert len(list(reader)) > 0
        
        with open(sites_file, 'r') as f:
            reader = csv.DictReader(f)
            assert len(list(reader)) > 0
