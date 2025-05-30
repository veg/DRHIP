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
            assert 'busted_pvalue' in row
            assert 'busted_lrt' in row
            assert 'busted_omega3' in row
            assert 'busted_evidence' in row
            
            # FEL results
            assert 'fel_sites_tested' in row
            assert 'fel_sites_positive_selection' in row
            assert 'fel_sites_negative_selection' in row
            
            # MEME results
            assert 'meme_sites_tested' in row
            assert 'meme_sites_selection' in row
            
            # Check data types
            assert isinstance(float(row['busted_pvalue']), float)
            assert isinstance(float(row['busted_lrt']), float)
            assert isinstance(float(row['busted_omega3']), float)

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
            
            # Check FEL data
            assert 'fel_alpha' in first_site
            assert 'fel_beta' in first_site
            assert 'fel_pvalue' in first_site
            assert 'fel_selection' in first_site
            
            # Check MEME data
            assert 'meme_alpha' in first_site
            assert 'meme_beta_neg' in first_site
            assert 'meme_beta_plus' in first_site
            assert 'meme_pvalue' in first_site
            
            # Check data types
            for field in ['fel_alpha', 'fel_beta', 'fel_pvalue']:
                assert isinstance(float(first_site[field]), float)

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
