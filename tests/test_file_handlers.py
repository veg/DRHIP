"""
Tests for file handling utilities.
"""

import os
import tempfile
import json
import pytest
from unittest.mock import patch, mock_open

from drhip.utils import file_handlers as fh
from drhip.methods.registry import HyPhyMethodRegistry


def test_get_genes():
    """Test get_genes function with mock directory structure."""
    # Create a temporary directory structure with mock results
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create method directories
        busted_dir = os.path.join(temp_dir, "BUSTED")
        fel_dir = os.path.join(temp_dir, "FEL")
        os.makedirs(busted_dir)
        os.makedirs(fel_dir)

        # Create mock result files
        with open(os.path.join(busted_dir, "gene1.BUSTED.json"), "w") as f:
            f.write("{}")
        with open(os.path.join(busted_dir, "gene2.BUSTED.json"), "w") as f:
            f.write("{}")
        with open(os.path.join(fel_dir, "gene1.FEL.json"), "w") as f:
            f.write("{}")
        with open(os.path.join(fel_dir, "gene3.FEL.json"), "w") as f:
            f.write("{}")

        # Test get_genes function
        genes = fh.get_genes(temp_dir)

        # Verify results
        assert isinstance(genes, list)
        assert set(genes) == {"gene1", "gene2", "gene3"}
        assert genes == sorted(genes)  # Check that genes are sorted


def test_get_genes_empty_directory():
    """Test get_genes function with an empty directory."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Test with empty directory
        with pytest.raises(FileNotFoundError):
            fh.get_genes(temp_dir)


def test_load_json_valid():
    """Test loading a valid JSON file."""
    # Create a temporary JSON file
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        json.dump({"key": "value"}, temp_file)
        temp_path = temp_file.name

    try:
        # Test loading the file
        result = fh.load_json(temp_path)

        # Verify result
        assert result == {"key": "value"}
    finally:
        # Clean up
        os.unlink(temp_path)


def test_load_json_nonexistent():
    """Test loading a nonexistent JSON file."""
    result = fh.load_json("/nonexistent/file.json")
    assert result is None


def test_load_json_invalid():
    """Test loading an invalid JSON file."""
    # Create a temporary file with invalid JSON
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        temp_file.write("{invalid: json")
        temp_path = temp_file.name

    try:
        # Test loading the file
        result = fh.load_json(temp_path)

        # Verify result
        assert result is None
    finally:
        # Clean up
        os.unlink(temp_path)
