"""
PRIME (PRoperty Informed Models of Evolution) method implementation.
"""

from typing import Dict, Any, List

from .base import HyPhyMethod

class PrimeMethod(HyPhyMethod):
    """Implementation of PRIME analysis processing."""
    
    PROPERTIES = [
        'Chemical composition',
        'Coil propensity',
        'Compressibility',
        'Hydropathy',
        'Molecular volume',
        'Polarity',
        'Refractivity',
        'Surface composition'
    ]
    
    def __init__(self):
        """Initialize PRIME method."""
        super().__init__("PRIME", "PRIME.json")
    
    def validate_input_json(self, results: Dict[str, Any]) -> List[str]:
        missing = []
        missing.extend(self.validate_required_paths(results, ['MLE.headers', 'MLE.content.0']))
        # If headers exist, ensure key columns are present
        missing_headers = self.missing_mle_headers(results, ['p-value'])
        if missing_headers:
            missing.extend([f"MLE.headers:{h}" for h in missing_headers])
        return missing
    
    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process PRIME results.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Processed results with standardized keys
        """
        processed = {}
        
        # Process sites under selection for each property
        for prop in self.PROPERTIES:
            sites_conserved = 0
            sites_altered = 0
            
            if self.has_mle_content(results) and self.has_mle_headers(results):
                # Get header indices
                header_indices = self.get_header_indices(results)
                
                # Get index for p-value
                pvalue_index = self.get_column_index(header_indices, 'p-value', 9)
                
                for row in results['MLE']['content']['0']:
                    p_value = float(row[pvalue_index])
                    if p_value <= 0.05:
                        sites_conserved += 1
                    else:
                        sites_altered += 1
            
            # We don't need to store property-specific counts in the summary
        
        return processed
    
    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific PRIME data.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Dictionary with site-specific metrics
        """
        # First, extract property tags and p-value indices
        prime_tags = []
        p_indices = []
        
        if self.has_mle_headers(results):
            headers = results['MLE']['headers']
            
            for i, (short_label, description) in enumerate(headers):
                if short_label == 'p-value':
                    p_indices.append(i)
                    prime_tags.append((0, 'overall'))
                elif short_label.startswith('p') and short_label[1:].isdigit():
                    p_indices.append(i)
                    # Extract property name from description
                    property_name = description.split(' ')[-1] if ' ' in description else short_label
                    prime_tags.append((int(short_label[1:]), property_name))
            
            # Sort tags by their index
            prime_tags = [tag for _, tag in sorted(prime_tags, key=lambda x: x[0])]
        
        # Define column names mapping - we'll handle p-value columns separately
        column_names = {}
        
        # Define a function to process each row
        def process_row(site_idx, row, column_indices):
            # Get p-values for all properties, with error handling
            pvals = []
            has_valid_data = True
            
            for idx in p_indices:
                try:
                    if idx < len(row):
                        pvals.append(float(row[idx]))
                    else:
                        has_valid_data = False  # Mark as invalid if index out of range
                        break
                except (ValueError, TypeError):
                    has_valid_data = False  # Mark as invalid if conversion fails
                    break
            
            # Determine the marker value
            if not has_valid_data or not pvals or not prime_tags:
                # For missing or malformed data
                return {
                    'prime_marker': 'NA'
                }
            elif min(pvals) <= 0.05:
                # For sites with significant properties
                significant_tags = [prime_tags[i] for i, pval in enumerate(pvals) if pval <= 0.05]
                return {
                    'prime_marker': ','.join(significant_tags)
                }
            else:
                # For non-significant sites, use dash as in original End2End pipeline
                return {
                    'prime_marker': '-'
                }
        
        # Use the helper to process site data
        return self.process_site_mle_data(results, column_names, process_row)
    
    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return []  # No summary fields needed
    
    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return [
            'prime_marker'
        ]
        
    @staticmethod
    def get_comparison_group_fields() -> List[str]:
        """Get list of fields that are specific to comparison groups."""
        return []  # PRIME doesn't have comparison group-specific fields
        
    def process_comparison_site_data(self, results: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from PRIME results.
        
        Args:
            results: Raw PRIME results dictionary
            
        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # PRIME doesn't have comparison group-specific data
        return {}
