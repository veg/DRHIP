"""
MEME (Mixed Effects Model of Evolution) method implementation.
"""

from typing import Any, Dict, List

from .base import HyPhyMethod


class MemeMethod(HyPhyMethod):
    """Implementation of MEME analysis processing."""

    def __init__(self):
        """Initialize MEME method."""
        super().__init__("MEME", "MEME.json")

    def validate_input_json(self, results: Dict[str, Any]) -> List[str]:
        missing = []
        missing.extend(
            self.validate_required_paths(results, ["MLE.headers", "MLE.content.0"])
        )
        missing_headers = self.missing_mle_headers(results, ["p-value"])
        if missing_headers:
            missing.extend([f"MLE.headers:{h}" for h in missing_headers])
        return missing

    def process_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process MEME results.

        Args:
            results: Raw MEME results dictionary

        Returns:
            Processed results with standardized keys
        """
        processed = {}

        # Process sites under selection
        positive_sites = 0
        if self.has_mle_content(results) and self.has_mle_headers(results):
            # Get header indices
            header_indices = self.get_header_indices(results)

            # Get index for p-value
            pvalue_index = self.get_column_index(
                header_indices, "p-value", 7
            )  # P-value for episodic selection

            for row in results["MLE"]["content"]["0"]:
                try:
                    p_value = float(row[pvalue_index])
                    if p_value <= 0.05:
                        positive_sites += 1
                except (ValueError, TypeError, IndexError):
                    # Skip rows with invalid data
                    pass

        # Store the count of sites under episodic selection in the summary
        processed["positive_sites"] = positive_sites

        return processed

    def process_site_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Process site-specific MEME data.

        Args:
            results: Raw MEME results dictionary

        Returns:
            Dictionary with site-specific metrics
        """
        # Define column names mapping
        column_names = {"site": "Site", "p-value": "p-value"}

        # Define a function to process each row
        def process_row(site_idx, row, column_indices):
            pvalue_index = column_indices["p-value"]

            # Check if we have valid indices and data
            if pvalue_index < 0 or pvalue_index >= len(row):
                return {"meme_marker": "NA"}  # Use NA for missing data

            try:
                pvalue = float(row[pvalue_index])

                # Set the marker based on significance
                if pvalue <= 0.05:
                    marker = f"{pvalue:.3f}"  # Format p-value for significant sites
                else:
                    marker = "-"  # Use dash for non-significant sites

                return {"meme_marker": marker}
            except (ValueError, TypeError):
                # Return NA for malformed data
                return {"meme_marker": "NA"}

        # Use the helper to process site data
        return self.process_site_mle_data(results, column_names, process_row)

    @staticmethod
    def get_summary_fields() -> List[str]:
        """Get list of summary fields produced by this method."""
        return ["positive_sites"]  # Return count of sites under episodic selection

    @staticmethod
    def get_site_fields() -> List[str]:
        """Get list of site-specific fields produced by this method."""
        return ["meme_marker"]

    @staticmethod
    def get_comparison_group_fields() -> List[str]:
        """Get list of fields that are specific to comparison groups."""
        return []  # MEME doesn't have comparison group-specific fields

    def process_comparison_site_data(
        self, results: Dict[str, Any]
    ) -> Dict[str, Dict[str, Any]]:
        """Process comparison group-specific site data from MEME results.

        Args:
            results: Raw MEME results dictionary

        Returns:
            Dictionary mapping site IDs to dictionaries of comparison group-specific data
        """
        # MEME doesn't have comparison group-specific data
        return {}
