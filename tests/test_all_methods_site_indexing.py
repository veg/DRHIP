"""
Tests for site indexing consistency across all HyPhy methods.

This test suite ensures that all methods use consistent site indexing (1-based)
and that site data is properly aligned when combining results from different methods.
"""

import os
import json
import unittest

from drhip.methods import HyPhyMethodRegistry
from drhip.parsers.process_gene import process_gene


class TestAllMethodsSiteIndexing(unittest.TestCase):
    """Test suite for site indexing consistency across all methods."""

    def setUp(self):
        """Set up test fixtures."""
        # Define paths to test data
        self.test_data_dir = os.path.join(os.path.dirname(__file__), 'data', 'hyphy')
        self.output_dir = os.path.join(os.path.dirname(__file__), 'output')
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize method instances using registry
        registry = HyPhyMethodRegistry()
        self.method_instances = {}
        
        # Create a mapping of method names to instances
        for method in registry.get_all_methods():
            method_name = method.name.upper()
            self.method_instances[method_name] = method
        
        # Map of method names to file paths
        self.method_files = {}
        
        # Find available test files for each method
        registry = HyPhyMethodRegistry()
        for method in registry.get_all_methods():
            method_name = method.name.upper()
            method_dir = os.path.join(self.test_data_dir, method_name)
            if os.path.exists(method_dir):
                for file in os.listdir(method_dir):
                    if file.endswith(f"{method_name}.json"):
                        self.method_files[method_name] = os.path.join(method_dir, file)
    
    def test_method_site_indexing(self):
        """Test that all methods use 1-based site indexing."""
        # Test each method that has site data
        registry = HyPhyMethodRegistry()
        
        # Only test methods that have site data (methods that implement get_site_fields)
        site_methods = [method.name.upper() for method in registry.get_all_methods() 
                       if hasattr(method, 'get_site_fields') and method.get_site_fields()]
        
        for method_name in site_methods:
            method_instance = self.method_instances.get(method_name)
            # Skip if no test file available
            if method_name not in self.method_files:
                continue
                
            # Load method results
            with open(self.method_files[method_name], 'r') as f:
                results = json.load(f)
            
            # Process site data
            site_data = method_instance.process_site_data(results)
            
            # Skip if method doesn't have site data
            if not site_data:
                continue
                
            # Check that site indices are 1-based
            site_indices = list(site_data.keys())
            
            # Verify we have site data
            self.assertTrue(
                len(site_indices) > 0, 
                f"No site data found in {method_name} results"
            )
            
            # Check that all site indices are integers and 1-based
            for site_idx in site_indices:
                self.assertIsInstance(
                    site_idx, 
                    int, 
                    f"Site index {site_idx} in {method_name} is not an integer"
                )
                self.assertGreaterEqual(
                    site_idx, 
                    1, 
                    f"Site index {site_idx} in {method_name} is not 1-based"
                )
            
            # Check that site 0 is not present (should be 1-based)
            self.assertNotIn(
                0, 
                site_data, 
                f"Site 0 found in {method_name} results (should be 1-based)"
            )
    
    def test_cfel_comparison_site_indexing(self):
        """Test that CFEL comparison site data uses 1-based indexing."""
        # Skip if no CFEL test file available
        if 'CFEL' not in self.method_files:
            return
            
        # Load CFEL results
        with open(self.method_files['CFEL'], 'r') as f:
            results = json.load(f)
        
        # Process comparison site data
        comparison_site_data = self.cfel.process_comparison_site_data(results)
        
        # Skip if no comparison site data
        if not comparison_site_data:
            return
            
        # Check each comparison group
        for group, site_data in comparison_site_data.items():
            # Check that site indices are 1-based
            site_indices = list(site_data.keys())
            
            # Verify we have site data
            self.assertTrue(
                len(site_indices) > 0, 
                f"No site data found in CFEL comparison group {group}"
            )
            
            # Check that all site indices are integers and 1-based
            for site_idx in site_indices:
                self.assertIsInstance(
                    site_idx, 
                    int, 
                    f"Site index {site_idx} in CFEL comparison group {group} is not an integer"
                )
                self.assertGreaterEqual(
                    site_idx, 
                    1, 
                    f"Site index {site_idx} in CFEL comparison group {group} is not 1-based"
                )
            
            # Check that site 0 is not present (should be 1-based)
            self.assertNotIn(
                0, 
                site_data, 
                f"Site 0 found in CFEL comparison group {group} (should be 1-based)"
            )
    
    def test_site_count_consistency(self):
        """Test that all methods have the same number of sites."""
        # Get site data from each method
        site_data_by_method = {}
        
        registry = HyPhyMethodRegistry()
        
        # Only test methods that have site data (methods that implement get_site_fields)
        site_methods = [method.name.upper() for method in registry.get_all_methods() 
                       if hasattr(method, 'get_site_fields') and method.get_site_fields()]
        
        for method_name in site_methods:
            method_instance = self.method_instances.get(method_name)
            # Skip if no test file available
            if method_name not in self.method_files:
                continue
                
            # Load method results
            with open(self.method_files[method_name], 'r') as f:
                results = json.load(f)
            
            # Process site data
            site_data = method_instance.process_site_data(results)
            
            # Skip if method doesn't have site data
            if not site_data:
                continue
                
            site_data_by_method[method_name] = site_data
        
        # Skip if fewer than 2 methods have site data
        if len(site_data_by_method) < 2:
            return
            
        # Check that all methods have the same number of sites
        site_counts = {method: len(site_data) for method, site_data in site_data_by_method.items()}
        reference_method = list(site_counts.keys())[0]
        reference_count = site_counts[reference_method]
        
        for method, count in site_counts.items():
            self.assertEqual(
                count, 
                reference_count,
                f"{method} has {count} sites, but {reference_method} has {reference_count} sites"
            )
    
    def test_site_index_alignment(self):
        """Test that site indices align across all methods."""
        # Get site data from each method
        site_data_by_method = {}
        
        registry = HyPhyMethodRegistry()
        
        # Only test methods that have site data (methods that implement get_site_fields)
        site_methods = [method.name.upper() for method in registry.get_all_methods() 
                       if hasattr(method, 'get_site_fields') and method.get_site_fields()]
        
        for method_name in site_methods:
            method_instance = self.method_instances.get(method_name)
            # Skip if no test file available
            if method_name not in self.method_files:
                continue
                
            # Load method results
            with open(self.method_files[method_name], 'r') as f:
                results = json.load(f)
            
            # Process site data
            site_data = method_instance.process_site_data(results)
            
            # Skip if method doesn't have site data
            if not site_data:
                continue
                
            site_data_by_method[method_name] = site_data
        
        # Skip if fewer than 2 methods have site data
        if len(site_data_by_method) < 2:
            return
            
        # Check that all methods have the same site indices
        reference_method = list(site_data_by_method.keys())[0]
        reference_indices = set(site_data_by_method[reference_method].keys())
        
        for method, site_data in site_data_by_method.items():
            if method == reference_method:
                continue
                
            method_indices = set(site_data.keys())
            
            # Check for missing indices in either method
            missing_in_reference = method_indices - reference_indices
            missing_in_method = reference_indices - method_indices
            
            self.assertEqual(
                len(missing_in_reference), 
                0, 
                f"Sites {missing_in_reference} are in {method} but missing in {reference_method}"
            )
            
            self.assertEqual(
                len(missing_in_method), 
                0, 
                f"Sites {missing_in_method} are in {reference_method} but missing in {method}"
            )
    
    def test_process_gene_site_indexing(self):
        """Test that process_gene produces correctly indexed site data."""
        # Skip if no test files available
        if not self.method_files:
            return
            
        # Get the gene name from the first available file
        gene_name = None
        for method_name, file_path in self.method_files.items():
            file_name = os.path.basename(file_path)
            gene_name = file_name.split(f".{method_name}")[0]
            break
            
        if not gene_name:
            return
            
        # Process gene with actual data
        process_gene(gene_name, self.test_data_dir, self.output_dir)
        
        # Check that the combined sites file was created
        combined_sites_file = os.path.join(self.output_dir, f'{gene_name}_sites.csv')
        self.assertTrue(
            os.path.exists(combined_sites_file), 
            f"Combined sites file not created: {combined_sites_file}"
        )
        
        # Read the combined sites file
        import pandas as pd
        sites_df = pd.read_csv(combined_sites_file)
        
        # Check that site indices are 1-based
        site_indices = sites_df['site'].tolist()
        
        # Check that site indices are sequential and 1-based
        self.assertEqual(min(site_indices), 1, "Site indices do not start at 1")
        
        # Check that there are no gaps in site indices
        expected_indices = set(range(1, max(site_indices) + 1))
        actual_indices = set(site_indices)
        
        self.assertEqual(
            expected_indices, 
            actual_indices, 
            f"Missing site indices: {expected_indices - actual_indices}"
        )
    
    def tearDown(self):
        """Clean up test fixtures."""
        # Remove output files
        import shutil
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)


if __name__ == '__main__':
    unittest.main()
