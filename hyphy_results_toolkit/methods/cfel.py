"""
CFEL (Contrast-FEL) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class CfelMethod(HyPhyMethod):
    """Implementation of CFEL analysis processing."""
    
    def __init__(self):
        """Initialize CFEL method."""
        # CFEL files are named as gene.CFEL.json in the contrastFEL directory
        super().__init__("CFEL", "CFEL.json")
        self._header_map = {}
        self._beta_idx_map = {}
        self._subs_idx_map = {}
    
    def _build_column_maps(self, results: Dict[str, Any], comparison_groups: List[str]) -> None:
        """Build column lookup maps for the CFEL data using base helper functions.
        
        Args:
            results: Raw results dictionary
            comparison_groups: List of comparison group names
        """
        # Use base helper to get header indices
        header_indices = self.get_header_indices(results)
        
        # Store the header map for use in other methods
        self._header_map = header_indices
        
        # Build maps for beta and substitution columns
        self._beta_idx_map = {g: self.get_column_index(header_indices, f"beta ({g})", -1) for g in comparison_groups}
        self._subs_idx_map = {g: self.get_column_index(header_indices, f"subs ({g})", -1) for g in comparison_groups}
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process CFEL results.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        # Initialize with default NA values for required fields
        processed = {
            'nt_conserved': 'NA',
            'aa_conserved': 'NA'
        }
        
        # If results are missing or invalid, return the defaults
        if not results or not isinstance(results, dict):
            print(f"CFEL.py: Invalid or missing CFEL results for {self._comparison_groups}")
            return processed
        
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
        
        # Build column lookup maps
        headers = results["MLE"]["headers"]
        # Always use the comparison groups we determined
        self._build_column_maps(headers, list(by_type.keys()))
        
        # Calculate number of sites and identify conserved sites
        try:
            if "MLE" in results and "content" in results["MLE"] and "0" in results["MLE"]["content"]:
                data_rows = results["MLE"]["content"]["0"]
                
                # Count conserved sites (both nt and aa)
                nt_conserved = 0
                aa_conserved = 0
                
                # Build column lookup maps if not already built
                if not self._beta_idx_map:
                    comparison_groups = list(by_type.keys())
                    self._build_column_maps(results, comparison_groups)
                
                for row in data_rows:
                    # Get beta values for all groups
                    beta_values = [row[idx] for idx in self._beta_idx_map.values() if idx >= 0]
                    
                    # Nucleotide conservation: all beta values are 0
                    if all(beta == 0.0 for beta in beta_values) if beta_values else False:
                        nt_conserved += 1
                    
                    # Amino acid conservation: synonymous substitutions only
                    # This requires all beta values to be 0 but some substitutions to exist
                    subs_values = [row[idx] for idx in self._subs_idx_map.values() if idx >= 0]
                    has_subs = any(subs > 0 for subs in subs_values) if subs_values else False
                    
                    if all(beta == 0.0 for beta in beta_values) if beta_values else False and has_subs:
                        aa_conserved += 1
                
                processed['nt_conserved'] = nt_conserved
                processed['aa_conserved'] = aa_conserved
            else:
                processed['sites'] = 'NA'
        except (KeyError, TypeError, IndexError):
            processed['sites'] = 'NA'
        
        # For each comparison group, calculate N and T
        for group, branches in by_type.items():
            # Number of sequences in this group
            processed[f'N_{group}'] = len(branches)
            
            # Total branch length for this group
            try:
                branch_lengths = [results['branch attributes']['0'][bn]['Global MG94xREV'] for bn in branches]
                processed[f'T_{group}'] = sum(branch_lengths)
            except (KeyError, TypeError):
                processed[f'T_{group}'] = 0.0
        
        # If there's just one group, use its values for the standard fields
        if len(by_type) == 1:
            group = list(by_type.keys())[0]
            processed['N'] = processed[f'N_{group}']
            processed['T'] = processed[f'T_{group}']
        
        # Get dN/dS for each group
        for k, r in results['fits']['Global MG94xREV']['Rate Distributions'].items():
            for group in by_type.keys():
                if group in k:
                    processed[f'dN/dS_{group}'] = r[0][0]
        
        # If there's just one group, use its dN/dS for the standard field
        if len(by_type) == 1:
            group = list(by_type.keys())[0]
            processed['dN/dS'] = processed[f'dN/dS_{group}']
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific CFEL data using base helper functions.
        
        Args:
            results: Raw CFEL results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        # Check if required data is available using base helper functions
        if not self.has_mle_content(results) or not self.has_mle_headers(results) or not results.get('tested', {}).get('0'):
            # If critical data is missing, return empty results
            return {}
        
        # Get the tag for each branch
        tested = results['tested']['0']
        by_type = {}
        
        # Require comparison groups to be set from process_gene.py
        if not self._comparison_groups:
            raise ValueError("Comparison groups must be set before processing CFEL site data")
            
        # Initialize by_type with empty lists for each comparison group
        for group in self._comparison_groups:
            by_type[group] = []
        
        # Now assign branches to the appropriate groups
        for branch, tag in tested.items():
            # Only add branches to groups that match our comparison groups
            if tag in by_type:
                by_type[tag].append(branch)
        
        # Build column lookup maps using base helper functions
        self._build_column_maps(results, list(by_type.keys()))
        
        # Define column name mapping for process_site_mle_data
        column_names = {
            'site': 'Site',
            'p-value': 'P-value (overall)'
        }
        
        # Get pairwise p-value column indices
        pairwise_pval_idx = {}
        for i, group1 in enumerate(self._comparison_groups):
            for group2 in self._comparison_groups[i+1:]:
                # Two ways the column might appear
                option_1 = f"P-value for {group1} vs {group2}"
                option_2 = f"P-value for {group2} vs {group1}"
                
                if option_1 in self._header_map:
                    pairwise_pval_idx[(group1, group2)] = self._header_map[option_1]
                elif option_2 in self._header_map:
                    pairwise_pval_idx[(group1, group2)] = self._header_map[option_2]
        
        # Define a function to process each row
        def process_row(site_idx, row, column_indices):
            site_data = {}
            
            # Get p-value index, with fallback
            p_overall_idx = column_indices.get('p-value', -1)
            if p_overall_idx == -1:
                # Try alternate name
                p_overall_idx = self._header_map.get("p-value", 5)  # Default to index 5 if not found
            
            # Check overall significance
            try:
                overall_pval = float(row[p_overall_idx])
                has_overall_significance = overall_pval <= 0.05
            except (ValueError, TypeError, IndexError):
                overall_pval = None
                has_overall_significance = False
            
            # Format marker similar to End2End-DENV implementation
            if overall_pval is None:
                marker = "NA"  # Use NA for missing or malformed data
            else:
                comparisons = []
                
                # Add overall p-value if significant
                if has_overall_significance:
                    comparisons.append(f"overall: {overall_pval:.3f}")
                    
                    # Check pairwise comparisons if overall is significant
                    for (group1, group2), idx in pairwise_pval_idx.items():
                        try:
                            pairwise_pval = float(row[idx])
                            if pairwise_pval <= 0.05:
                                # Format comparison key with groups in alphabetical order for consistency
                                sorted_groups = sorted([group1, group2])
                                comparisons.append(f"{sorted_groups[0]} vs {sorted_groups[1]}: {pairwise_pval:.3f}")
                        except (ValueError, TypeError, IndexError):
                            pass  # Skip this comparison if we can't get a valid p-value
                
                # Format marker with all significant comparisons
                marker = ", ".join(comparisons) if comparisons else "-"  # NS = Not Significant
            
            # Store the marker for this site
            site_data["marker"] = marker
            
            # Get beta values for each group
            for group, idx in self._beta_idx_map.items():
                try:
                    site_data[f"beta_{group}"] = float(row[idx]) if idx >= 0 else None
                except (ValueError, TypeError, IndexError):
                    site_data[f"beta_{group}"] = None
            
            # Get substitution counts for each group
            for group, idx in self._subs_idx_map.items():
                try:
                    site_data[f"subs_{group}"] = float(row[idx]) if idx >= 0 else None
                except (ValueError, TypeError, IndexError):
                    site_data[f"subs_{group}"] = None
            
            # Check if site is conserved (nucleotide or amino acid)
            try:
                # Get beta values for all groups
                beta_values = [float(row[idx]) for idx in self._beta_idx_map.values() if idx >= 0]
                
                # Nucleotide conservation: all beta values are 0
                is_nt_conserved = all(beta == 0.0 for beta in beta_values) if beta_values else False
                
                # Amino acid conservation: synonymous substitutions only (beta = 0 but alpha might be > 0)
                subs_values = [float(row[idx]) for idx in self._subs_idx_map.values() if idx >= 0]
                has_subs = any(subs > 0 for subs in subs_values) if subs_values else False
                is_aa_conserved = is_nt_conserved and has_subs
                
                site_data["nt_conserved"] = is_nt_conserved
                site_data["aa_conserved"] = is_aa_conserved
            except (IndexError, TypeError, ValueError):
                site_data["nt_conserved"] = False
                site_data["aa_conserved"] = False
            
            return site_data
        
        # Use the base helper to process site data
        return self.process_site_mle_data(results, column_names, process_row)

    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return [
            'nt_conserved',
            'aa_conserved'
        ]
        
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'cfel_marker',
            'nt_conserved',
            'aa_conserved'
        ]