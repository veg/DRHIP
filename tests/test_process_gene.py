"""
Tests for process_gene module.
"""

import os
import csv
import tempfile
import shutil
from hyphy_results_toolkit.parsers.process_gene import process_gene

def test_process_gene_output_files(results_dir):
    """Test that process_gene creates expected output files."""
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Process gene
        gene_name = 'capsid_protein_C'
        process_gene(gene_name, results_dir, output_dir)
        
        # Check that output files exist
        assert os.path.exists(os.path.join(output_dir, f'{gene_name}_summary.csv'))
        assert os.path.exists(os.path.join(output_dir, f'{gene_name}_sites.csv'))

def test_process_gene_summary_content(results_dir):
    """Test that process_gene produces correct summary content."""
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Process gene
        gene_name = 'capsid_protein_C'
        process_gene(gene_name, results_dir, output_dir)
        
        # Read summary file
        summary_file = os.path.join(output_dir, f'{gene_name}_summary.csv')
        with open(summary_file, 'r') as f:
            reader = csv.DictReader(f)
            row = next(reader)
            
            # Check that key fields exist and have valid values
            # BUSTED results
            assert 'BUSTED_pval' in row
            assert float(row['BUSTED_pval']) >= 0
            assert 'BUSTED_omega3' in row
            assert 'BUSTED_prop_sites_in_omega3' in row
            
            # FEL/MEME results
            assert 'N' in row
            assert int(row['N']) > 0
            assert 'positive_sites' in row
            assert 'negative_sites' in row

def test_process_gene_sites_content(results_dir):
    """Test that process_gene produces correct site-specific content."""
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Process gene
        gene_name = 'capsid_protein_C'
        process_gene(gene_name, results_dir, output_dir)
        
        # Read sites file
        sites_file = os.path.join(output_dir, f'{gene_name}_sites.csv')
        with open(sites_file, 'r') as f:
            reader = csv.DictReader(f)
            sites = list(reader)
            
            # Check that we have site data
            assert len(sites) > 0
            
            # Check that key fields exist in the first site
            first_site = sites[0]
            assert 'gene' in first_site
            assert 'site' in first_site
            
            # Check for some site-specific data fields
            # Field names may vary based on the current implementation
            # Just check for a few key fields that should be present
            assert 'gene' in first_site
            assert 'site' in first_site
            assert 'majority_residue' in first_site
            
            # Check for at least one method-specific field
            method_fields = [key for key in first_site.keys() 
                            if any(method in key.lower() for method in ['fel', 'meme', 'busted'])]                
            assert len(method_fields) > 0

def test_process_gene_comparison_site_data(results_dir):
    """Test that process_gene produces correct comparison group site-specific content with sequence analysis."""
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Process gene
        gene_name = 'capsid_protein_C'
        process_gene(gene_name, results_dir, output_dir)
        
        # Check if comparison site file exists
        comparison_sites_file = os.path.join(output_dir, f'{gene_name}_comparison_sites.csv')
        
        # If the file exists, check its contents
        if os.path.exists(comparison_sites_file):
            with open(comparison_sites_file, 'r') as f:
                reader = csv.DictReader(f)
                comp_sites = list(reader)
                
                # If we have comparison site data, check for sequence analysis fields
                if comp_sites:
                    first_comp_site = comp_sites[0]
                    
                    # Check for required fields
                    assert 'gene' in first_comp_site
                    assert 'site' in first_comp_site
                    assert 'comparison_group' in first_comp_site
                    
                    # Check for sequence analysis fields
                    # These may not all be present depending on the data
                    sequence_fields = ['unique_aas', 'has_diff_majority', 'majority_residue', 
                                      'aa_diversity', 'composition']
                    
                    # At least some of these fields should be present
                    present_seq_fields = [field for field in sequence_fields if field in first_comp_site]
                    assert len(present_seq_fields) > 0, "No sequence analysis fields found in comparison site data"

def test_process_gene_thread_safety(results_dir):
    """Test that process_gene is thread-safe when writing output."""
    import threading
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Create temporary directory with multiple copies of the same data
        # to simulate multiple genes
        with tempfile.TemporaryDirectory() as temp_data_dir:
            # Create subdirectories for each method
            for method in ['BUSTED', 'FEL', 'MEME', 'PRIME']:
                os.makedirs(os.path.join(temp_data_dir, method), exist_ok=True)
            
            # Copy the same data for multiple gene names
            genes = ['gene1', 'gene2', 'gene3']
            for gene in genes:
                for method in ['BUSTED', 'FEL', 'MEME', 'PRIME']:
                    source_dir = os.path.join(results_dir, method)
                    target_dir = os.path.join(temp_data_dir, method)
                    
                    # Get the first JSON file in the source directory
                    source_files = [f for f in os.listdir(source_dir) if f.endswith('.json')]
                    if source_files:
                        source_file = os.path.join(source_dir, source_files[0])
                        target_file = os.path.join(target_dir, f'{gene}.{method}.json')
                        shutil.copy(source_file, target_file)
            
            # Process genes in parallel
            threads = []
            for gene in genes:
                thread = threading.Thread(
                    target=process_gene,
                    args=(gene, temp_data_dir, output_dir)
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
