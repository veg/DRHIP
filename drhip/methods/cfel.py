"""
CFEL (Contrast-FEL) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class CfelMethod(HyPhyMethod):
    """Implementation of CFEL analysis processing."""
    
    def __init__(self):
        """Initialize CFEL method."""
        # CFEL files are named as gene.CONTRASTFEL.json in the CONTRASTFEL directory
        super().__init__("CFEL", "CONTRASTFEL.json")
        self._header_map = {}
        self._beta_idx_map = {}
        self._subs_idx_map = {}
    
    def _build_column_maps(self, results: Dict[str, Any], comparison_groups: List[str]) -> None:
        """Build column lookup maps for the CFEL data using base helper functions.
        
        Args:
            results: Raw results dictionary
            comparison_groups: List of comparison group names
        """
        # Skip if no MLE headers
        if not self.has_mle_headers(results):
            return
            
        # Use base helper to get header indices
        header_indices = self.get_header_indices(results)
        
        # Store the header map for use in other methods
        self._header_map = header_indices
        
        # Initialize maps
        self._beta_idx_map = {}
        self._subs_idx_map = {}
        
        # Find column indices for beta and substitution values
        for header_name, idx in header_indices.items():
            if 'beta' in header_name.lower():
                # Extract group name from column header
                import re
                match = re.search(r'beta\s*\(([^)]+)\)', header_name)
                if match:
                    group = match.group(1)
                    if group in comparison_groups:
                        self._beta_idx_map[group] = idx
            elif 'subs' in header_name.lower():
                # Extract group name from substitutions column
                match = re.search(r'subs\s*\(([^)]+)\)', header_name, re.IGNORECASE)
                if match:
                    group = match.group(1)
                    if group in comparison_groups:
                        self._subs_idx_map[group] = idx
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process CFEL results.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Dictionary with diff_sites count from CFEL analysis
        """
        # Initialize with empty dictionary
        processed = {}
        
        # Calculate diff_sites (sites with Q-value < 0.05 from CFEL MLE output)
        if self.has_mle_headers(results) and self.has_mle_content(results):
            try:
                # Get header indices
                header_indices = self.get_header_indices(results)
                
                # Find the Q-value column index
                q_value_idx = -1
                for header_name, idx in header_indices.items():
                    if 'q-value' in header_name.lower():
                        q_value_idx = idx
                        break
                
                # Count sites with Q-value < 0.05
                diff_sites_count = 0
                if q_value_idx >= 0:
                    for row in results['MLE']['content']['0']:
                        if q_value_idx < len(row):
                            try:
                                q_value = float(row[q_value_idx])
                                if q_value <= 0.05:  # Significant if Q-value <= 0.05
                                    diff_sites_count += 1
                            except (ValueError, TypeError):
                                # Skip rows with invalid Q-values
                                continue
                
                # Add diff_sites to processed results
                processed['diff_sites'] = diff_sites_count
                
            except Exception as e:
                print(f"Error calculating diff_sites from CFEL results: {e}")
        
        return processed
        
    def process_comparison_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process comparison group data that is not site-specific.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Dictionary with comparison group data
        """
        # If results are missing or invalid, return empty dict
        if not results or not isinstance(results, dict):
            print(f"CFEL.py: Invalid or missing CFEL results for {self._comparison_groups}")
            return {}
        
        # Get the tag for each branch
        tested = results['tested']['0']
        by_type: Dict[str, List[str]] = {}
        
        # Require comparison groups to be set from process_gene.py
        if not self._comparison_groups:
            raise ValueError("Comparison groups must be set before processing CFEL results")
            
        # Initialize by_type with empty lists for each comparison group
        for group in self._comparison_groups:
            by_type[group] = []
        
        # Now assign branches to the appropriate groups
        for branch, tag in tested.items():
            # Only add branches to groups that match our comparison groups
            if tag in by_type:
                by_type[tag].append(branch)
        
        # Build column lookup maps for site analysis
        if self.has_mle_headers(results):
            headers = results["MLE"]["headers"]
            # Always use the comparison groups we determined
            self._build_column_maps(headers, list(by_type.keys()))
        
        # Initialize comparison data
        comparison_data = {}
        
        # For each comparison group, calculate N and T
        for group, branches in by_type.items():
            group_data = {}
            
            # Number of sequences in this group
            group_data['group_N'] = len(branches)
            
            # Total branch length for this group
            try:
                branch_lengths = [results['branch attributes']['0'][bn]['Global MG94xREV'] for bn in branches]
                group_data['group_T'] = sum(branch_lengths)
            except (KeyError, TypeError):
                group_data['group_T'] = 0.0
            
            # Get dN/dS for this group if available
            try:
                for k, r in results['fits']['Global MG94xREV']['Rate Distributions'].items():
                    if group in k:
                        group_data['group_dN/dS'] = r[0][0]
                        break
            except (KeyError, TypeError, IndexError):
                group_data['group_dN/dS'] = 'NA'
            
            # Calculate conserved sites for this group
            if self.has_mle_content(results):
                try:
                    data_rows = results["MLE"]["content"]["0"]
                    
                    # Count conserved sites (both nt and aa)
                    nt_conserved = 0
                    aa_conserved = 0
                    
                    for row in data_rows:
                        # Get beta value for this group
                        beta_idx = self._beta_idx_map.get(group, -1)
                        if beta_idx >= 0:
                            beta = float(row[beta_idx])
                            
                            # Nucleotide conservation: beta = 0
                            is_nt_conserved = (beta == 0.0)
                            
                            # Amino acid conservation: beta = 0 but substitutions exist
                            subs_idx = self._subs_idx_map.get(group, -1)
                            has_subs = False
                            if subs_idx >= 0:
                                has_subs = (float(row[subs_idx]) > 0)
                            
                            is_aa_conserved = is_nt_conserved and has_subs
                            
                            if is_nt_conserved:
                                nt_conserved += 1
                            if is_aa_conserved:
                                aa_conserved += 1
                    
                    group_data['group_nt_conserved'] = nt_conserved
                    group_data['group_aa_conserved'] = aa_conserved
                except (KeyError, TypeError, IndexError, ValueError):
                    group_data['group_nt_conserved'] = 'NA'
                    group_data['group_aa_conserved'] = 'NA'
            
            comparison_data[group] = group_data
        
        return comparison_data
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific CFEL data using base helper functions.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Empty dictionary as CFEL doesn't use standard site fields
        """
        # CFEL doesn't use standard site fields, only comparison group site fields
        return {}
        
    def process_comparison_site_data(self, results: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from CFEL results.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # Skip processing if no MLE content or headers or no comparison groups
        if not self.has_mle_content(results) or not self.has_mle_headers(results) or not self._comparison_groups:
            print("Skipping CFEL comparison site data processing due to missing data")
            return {}
            
        # Get MLE content
        mle_content = results['MLE']['content']

        self._build_column_maps(results, list(self._comparison_groups))
        
        # Find column indices for beta values by group
        beta_idx_map = self._beta_idx_map
        
        # Get branch information for calculating N and T values
        tested = results.get('tested', {}).get('0', {})
        by_type = {group: [] for group in self._comparison_groups}
        
        # Assign branches to the appropriate groups
        for branch, tag in tested.items():
            if tag in by_type:
                by_type[tag].append(branch)
        
        # Calculate N and T for each comparison group
        group_N_values = {}
        group_T_values = {}
        
        for group, branches in by_type.items():
            # Number of sequences in this group
            group_N_values[group] = len(branches)
            
            # Total branch length for this group
            try:
                branch_lengths = [results['branch attributes']['0'][bn]['Global MG94xREV'] for bn in branches]
                group_T_values[group] = sum(branch_lengths)
            except (KeyError, TypeError):
                group_T_values[group] = 0.0
        
        # Initialize result dictionary
        comparison_data = {}
        
        # Process site data - in CONTRASTFEL format, '0' contains a list of rows
        if '0' in mle_content and isinstance(mle_content['0'], list):
            rows = mle_content['0']
            
            # Each row represents a site
            for site_idx, row in enumerate(rows):
                site_id = str(site_idx + 1)  # Convert to 1-based site index as string
                site_comparison_data = {}
                
                # Process data for each comparison group
                for group in self._comparison_groups:
                    group_data = {}
                    
                    # Add CFEL marker for this group
                    try:
                        beta_idx = beta_idx_map.get(group, -1)
                        if beta_idx >= 0 and beta_idx < len(row):
                            beta = float(row[beta_idx])
                            
                            # Set CFEL marker based on beta value
                            if beta > 0:
                                group_data['cfel_marker'] = '+'
                            elif beta < 0:
                                group_data['cfel_marker'] = '-'
                            else:
                                group_data['cfel_marker'] = '='
                        else:
                            group_data['cfel_marker'] = '?'
                    except (ValueError, TypeError, IndexError):
                        group_data['cfel_marker'] = '?'
                    
                    site_comparison_data[group] = group_data
                
                comparison_data[site_id] = site_comparison_data
        
        return comparison_data
        
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        # CFEL now provides the diff_sites field
        return ['diff_sites']
        
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        # CFEL doesn't use standard site fields
        return []
        
    @staticmethod
    def get_comparison_group_site_fields() -> List[str]:
        """Get list of site-specific fields that are specific to comparison groups."""
        return [
            'cfel_marker',  # CFEL marker for this site in this comparison group
        ]
        
    @staticmethod
    def get_comparison_group_summary_fields() -> List[str]:
        """Get list of non-site-specific fields that are specific to comparison groups."""
        return [
            'group_N',              # Number of sequences in this comparison group
            'group_T',              # Total branch length for this comparison group
            'group_dN/dS',          # dN/dS ratio for this comparison group
            'group_nt_conserved',   # Number of nucleotide conserved sites in this group
            'group_aa_conserved'    # Number of amino acid conserved sites in this group
        ]