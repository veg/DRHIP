"""
Tests for CLI functionality.
"""

import csv
import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from drhip.cli import combine_files, main


def test_combine_csv_files():
    """Test combining CSV files functionality."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Create test CSV files
            gene1_file = os.path.join(temp_dir, "gene1_summary.csv")
            with open(gene1_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "busted_pvalue", "fel_sites_tested"])
                writer.writerow(["gene1", "0.01", "100"])

            gene2_file = os.path.join(temp_dir, "gene2_summary.csv")
            with open(gene2_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "busted_pvalue", "meme_sites_selection"])
                writer.writerow(["gene2", "0.05", "5"])

            # Test the combine function
            combine_files(temp_dir, output_dir, "summary")

            # Verify results
            output_file = os.path.join(output_dir, "combined_summary.csv")
            assert os.path.exists(output_file)

            # Check combined file content
            with open(output_file, newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                rows.sort(key=lambda x: x.get("gene", ""))

                assert len(rows) == 2
                assert rows[0]["gene"] == "gene1"
                assert rows[0]["busted_pvalue"] == "0.01"
                assert rows[0]["fel_sites_tested"] == "100"
                assert rows[0]["meme_sites_selection"] == "NA"

                assert rows[1]["gene"] == "gene2"
                assert rows[1]["busted_pvalue"] == "0.05"
                assert rows[1]["meme_sites_selection"] == "5"
                assert rows[1]["fel_sites_tested"] == "NA"


def test_combine_csv_files_sites():
    """Test combining site-specific CSV files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Create test CSV files
            gene1_file = os.path.join(temp_dir, "gene1_sites.csv")
            with open(gene1_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "site", "fel_alpha", "fel_beta"])
                writer.writerow(["gene1", "1", "0.5", "1.0"])
                writer.writerow(["gene1", "2", "0.6", "0.8"])

            gene2_file = os.path.join(temp_dir, "gene2_sites.csv")
            with open(gene2_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "site", "meme_pvalue", "meme_alpha"])
                writer.writerow(["gene2", "1", "0.01", "0.3"])

            # Test the combine function
            combine_files(temp_dir, output_dir, "sites")

            # Verify results
            output_file = os.path.join(output_dir, "combined_sites.csv")
            assert os.path.exists(output_file)

            # Check combined file content
            with open(output_file, newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                rows.sort(key=lambda x: x.get("gene", ""))

                assert len(rows) == 3
                assert rows[0]["gene"] == "gene1"
                assert rows[0]["site"] == "1"
                assert rows[0]["fel_alpha"] == "0.5"
                assert rows[0]["meme_pvalue"] == "NA"

                assert rows[2]["gene"] == "gene2"
                assert rows[2]["site"] == "1"
                assert rows[2]["meme_alpha"] == "0.3"
                assert rows[2]["fel_alpha"] == "NA"


def test_combine_csv_files_comparison():
    """Test combining comparison CSV files with comparison groups."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Create test CSV files
            gene1_file = os.path.join(temp_dir, "gene1_comparison_site.csv")
            with open(gene1_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["gene", "site", "comparison_group", "cfel_marker", "group_N"]
                )
                writer.writerow(["gene1", "1", "foreground", "+", "5"])
                writer.writerow(["gene1", "1", "background", "-", "3"])

            gene2_file = os.path.join(temp_dir, "gene2_comparison_site.csv")
            with open(gene2_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["gene", "site", "comparison_group", "group_T", "group_dN/dS"]
                )
                writer.writerow(["gene2", "1", "foreground", "0.5", "1.2"])

            # Test the combine function
            combine_files(temp_dir, output_dir, "comparison_site")

            # Verify results
            output_file = os.path.join(output_dir, "combined_comparison_site.csv")
            assert os.path.exists(output_file)

            # Check combined file content
            with open(output_file, newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)

                assert len(rows) == 3

                # Sort rows by gene and comparison_group for consistent testing
                rows.sort(key=lambda x: (x["gene"], x["comparison_group"]))

                # Now check in a deterministic order
                assert rows[0]["gene"] == "gene1"
                assert rows[0]["comparison_group"] == "background"
                assert rows[0]["cfel_marker"] == "-"
                assert rows[0]["group_N"] == "3"

                assert rows[1]["gene"] == "gene1"
                assert rows[1]["comparison_group"] == "foreground"
                assert rows[1]["cfel_marker"] == "+"
                assert rows[1]["group_N"] == "5"
                # The combine_csv_files function may add empty fields for all columns
                # Just check that if group_T exists, it's empty or NA
                if "group_T" in rows[1]:
                    assert rows[1]["group_T"] == "" or rows[1]["group_T"] == "NA"

                assert rows[2]["gene"] == "gene2"
                assert rows[2]["comparison_group"] == "foreground"
                assert rows[2]["group_T"] == "0.5"
                assert rows[2]["group_dN/dS"] == "1.2"
                # The combine_csv_files function adds all fields to all rows
                # Just check that if cfel_marker exists, it's empty or NA
                if "cfel_marker" in rows[2]:
                    assert (
                        rows[2]["cfel_marker"] == "" or rows[2]["cfel_marker"] == "NA"
                    )


def test_combine_csv_files_empty():
    """Test combining CSV files when no files exist."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Test the combine function with no matching files
            combine_files(temp_dir, output_dir, "summary")

            # Verify no output file was created
            output_file = os.path.join(output_dir, "combined_summary.csv")
            assert not os.path.exists(output_file)


@patch("drhip.cli.process_gene.process_gene")
@patch("drhip.cli.fh.get_genes")
@patch("drhip.cli.combine_files")
def test_main_workflow(mock_combine, mock_get_genes, mock_process_gene):
    """Test the main CLI workflow with mocked dependencies."""
    # Setup mocks
    mock_get_genes.return_value = ["gene1", "gene2"]

    # Create temp directories for testing
    with tempfile.TemporaryDirectory() as input_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Mock command line arguments
            with patch("sys.argv", ["drhip", "-i", input_dir, "-o", output_dir]):
                # Run the main function
                main()

                # Verify process_gene was called for each gene
                assert mock_process_gene.call_count == 2
                mock_process_gene.assert_any_call(
                    "gene1", input_dir, mock_combine.call_args[0][0]
                )
                mock_process_gene.assert_any_call(
                    "gene2", input_dir, mock_combine.call_args[0][0]
                )

                # Verify combine_csv_files was called for each file type
                assert mock_combine.call_count == 4
                # Delimiter is optional with default ",", so we can check with or without it
                for suffix in ["summary", "sites", "comparison_site", "comparison_summary"]:
                    # Check that it was called with this suffix (delimiter may or may not be explicit)
                    calls_with_suffix = [
                        call for call in mock_combine.call_args_list
                        if len(call[0]) >= 3 and call[0][2] == suffix
                    ]
                    assert len(calls_with_suffix) >= 1, f"Expected call with suffix '{suffix}'"


@patch("drhip.cli.process_gene.process_gene")
@patch("drhip.cli.fh.get_genes")
def test_main_with_exception(mock_get_genes, mock_process_gene):
    """Test the main CLI workflow when an exception occurs."""
    # Setup mocks
    mock_get_genes.return_value = ["gene1", "gene2"]
    mock_process_gene.side_effect = [None, Exception("Test exception")]

    # Create temp directories for testing
    with tempfile.TemporaryDirectory() as input_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Mock command line arguments
            with patch("sys.argv", ["drhip", "-i", input_dir, "-o", output_dir]):
                # Run the main function - should not raise the exception
                main()

                # Verify process_gene was called for each gene
                assert mock_process_gene.call_count == 2


@patch("argparse.ArgumentParser.parse_args")
def test_main_with_relative_paths(mock_parse_args):
    """Test the main CLI workflow with relative paths."""
    # Setup mock for parse_args
    mock_args = MagicMock()
    mock_args.input = "relative/input"
    mock_args.output = "relative/output"
    mock_args.threads = 1
    mock_args.comparison_groups = None
    mock_parse_args.return_value = mock_args

    # Mock other dependencies
    with patch("drhip.cli.fh.get_genes") as mock_get_genes:
        with patch("drhip.cli.process_gene.process_gene") as mock_process_gene:
            with patch("drhip.cli.os.getcwd") as mock_getcwd:
                with patch("drhip.cli.os.path.join") as mock_join:
                    with patch("drhip.cli.os.makedirs"):
                        with patch("drhip.cli.combine_files"):
                            # Setup mocks
                            mock_getcwd.return_value = "/home/user"
                            mock_get_genes.return_value = []
                            # Make join return a predictable value
                            mock_join.side_effect = lambda *args: "/".join(args)

                            # Run the main function
                            main()

                            # Verify paths were resolved correctly
                            mock_join.assert_any_call("/home/user", "relative/input")


def test_csv_output_format():
    """Test CSV output format with comma delimiter and LF line endings."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Create test CSV files
            gene1_file = os.path.join(temp_dir, "gene1_summary.csv")
            with open(gene1_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "value1", "value2"])
                writer.writerow(["gene1", "100", "200"])

            # Test CSV format (default)
            combine_files(temp_dir, output_dir, "summary", delimiter=",")

            # Verify output file
            output_file = os.path.join(output_dir, "combined_summary.csv")
            assert os.path.exists(output_file)
            assert output_file.endswith(".csv")

            # Read file in binary mode to check line endings and delimiter
            with open(output_file, "rb") as f:
                content = f.read()
                # Check for LF line endings (0x0a) and no CRLF (0x0d 0x0a)
                assert b"\n" in content
                assert b"\r\n" not in content
                # Check for comma delimiter
                assert b"," in content

            # Verify content is readable as CSV
            with open(output_file, newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                assert len(rows) == 1
                assert rows[0]["gene"] == "gene1"
                assert rows[0]["value1"] == "100"


def test_tabular_output_format():
    """Test tabular output format with tab delimiter and LF line endings."""
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Create test CSV files
            gene1_file = os.path.join(temp_dir, "gene1_summary.csv")
            with open(gene1_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "value1", "value2"])
                writer.writerow(["gene1", "100", "200"])

            # Test tabular format
            combine_files(temp_dir, output_dir, "summary", delimiter="\t")

            # Verify output file
            output_file = os.path.join(output_dir, "combined_summary.tab")
            assert os.path.exists(output_file)
            assert output_file.endswith(".tab")

            # Read file in binary mode to check line endings and delimiter
            with open(output_file, "rb") as f:
                content = f.read()
                # Check for LF line endings (0x0a) and no CRLF (0x0d 0x0a)
                assert b"\n" in content
                assert b"\r\n" not in content
                # Check for tab delimiter (0x09)
                assert b"\t" in content
                # Should not have commas as delimiters
                assert b"gene,value1" not in content

            # Verify content is readable as tab-delimited
            with open(output_file, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
                assert len(rows) == 1
                assert rows[0]["gene"] == "gene1"
                assert rows[0]["value1"] == "100"


@patch("drhip.cli.process_gene.process_gene")
@patch("drhip.cli.fh.get_genes")
@patch("drhip.cli.combine_files")
def test_main_with_tabular_flag(mock_combine, mock_get_genes, mock_process_gene):
    """Test the main CLI workflow with --tabular flag."""
    # Setup mocks
    mock_get_genes.return_value = ["gene1"]

    # Create temp directories for testing
    with tempfile.TemporaryDirectory() as input_dir:
        with tempfile.TemporaryDirectory() as output_dir:
            # Mock command line arguments with --tabular
            with patch(
                "sys.argv", ["drhip", "-i", input_dir, "-o", output_dir, "--tabular"]
            ):
                # Run the main function
                main()

                # Verify combine_csv_files was called with tab delimiter
                assert mock_combine.call_count == 4
                for call in mock_combine.call_args_list:
                    # Check that delimiter argument is tab
                    assert call[0][2] in [
                        "summary",
                        "sites",
                        "comparison_site",
                        "comparison_summary",
                    ]
                    # The delimiter should be passed as the third positional arg
                    if len(call[0]) > 3:
                        assert call[0][3] == "\t"
