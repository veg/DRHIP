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

def test_process_gene_comparison_data(comparison_results_dir):
    """Test that process_gene produces correct comparison group data with the comparison test dataset."""
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as output_dir:
        # Process gene - use the full gene name as it appears in the comparison test data files
        gene_name = 'pretend_DENV1_ref.part_NC_001477.1__capsid_protein_C__95-394_DENV1'
        process_gene(gene_name, comparison_results_dir, output_dir)
        
        # Check that standard output files exist
        assert os.path.exists(os.path.join(output_dir, f'{gene_name}_summary.csv'))
        assert os.path.exists(os.path.join(output_dir, f'{gene_name}_sites.csv'))
        
        # Check that comparison-specific output files exist
        comparison_site_file = os.path.join(output_dir, f'{gene_name}_comparison_site.csv')
        comparison_summary_file = os.path.join(output_dir, f'{gene_name}_comparison_summary.csv')
        
        assert os.path.exists(comparison_site_file), "Comparison site file not created"
        assert os.path.exists(comparison_summary_file), "Comparison summary file not created"
        
        # Validate comparison site file contents
        with open(comparison_site_file, 'r') as f:
            reader = csv.DictReader(f)
            comp_sites = list(reader)
            
            # Check that we have site data
            assert len(comp_sites) > 0, "No comparison site data found"
            
            # Check the first site for required fields
            first_comp_site = comp_sites[0]
            assert 'gene' in first_comp_site
            assert 'site' in first_comp_site
            assert 'comparison_group' in first_comp_site
            
            # Check for method-specific fields from CONTRASTFEL or RELAX
            method_fields = [key for key in first_comp_site.keys() 
                           if any(method in key.lower() for method in ['cfel', 'relax'])]
            assert len(method_fields) > 0, "No method-specific fields found in comparison site data"
        
        # Validate comparison summary file contents
        with open(comparison_summary_file, 'r') as f:
            reader = csv.DictReader(f)
            comp_summaries = list(reader)
            
            # Check that we have summary data
            assert len(comp_summaries) > 0, "No comparison summary data found"
            
            # Check the first summary for required fields
            first_comp_summary = comp_summaries[0]
            assert 'gene' in first_comp_summary
            assert 'comparison_group' in first_comp_summary
            
            # Check for expected comparison group fields from CFEL
            expected_fields = ['group_N', 'group_T', 'group_dN/dS']
            found_fields = [field for field in expected_fields if field in first_comp_summary]
            print(f"Found comparison group fields: {found_fields}")
            assert len(found_fields) > 0, "No expected comparison group fields found in comparison summary data"

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
