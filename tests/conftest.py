"""
Pytest configuration and fixtures.
"""

import os
import json
from pathlib import Path
from typing import Dict, Any
import pytest

@pytest.fixture
def results_dir() -> str:
    """Get the path to the HyPhy results directory for tests."""
    return str(Path(__file__).parent / 'data' / 'hyphy')

@pytest.fixture
def comparison_results_dir() -> str:
    """Get the path to the HyPhy comparison results directory for tests."""
    return str(Path(__file__).parent / 'data' / 'hyphy-comparisons')

# Fixtures for real HyPhy test data
@pytest.fixture
def real_data_dir() -> str:
    """Get the path to the real HyPhy test data directory."""
    return str(Path(__file__).parent / 'data' / 'hyphy')

@pytest.fixture
def real_fel_results(real_data_dir) -> Dict[str, Any]:
    """Load real FEL results from test data."""
    file_path = os.path.join(real_data_dir, 'FEL', 'capsid_protein_C.FEL.json')
    with open(file_path, 'r') as f:
        return json.load(f)

@pytest.fixture
def real_busted_results(real_data_dir) -> Dict[str, Any]:
    """Load real BUSTED results from test data."""
    file_path = os.path.join(real_data_dir, 'BUSTED', 'capsid_protein_C.BUSTED.json')
    with open(file_path, 'r') as f:
        return json.load(f)

@pytest.fixture
def real_meme_results(real_data_dir) -> Dict[str, Any]:
    """Load real MEME results from test data."""
    file_path = os.path.join(real_data_dir, 'MEME', 'capsid_protein_C.MEME.json')
    with open(file_path, 'r') as f:
        return json.load(f)

@pytest.fixture
def real_prime_results(real_data_dir) -> Dict[str, Any]:
    """Load real PRIME results from test data."""
    file_path = os.path.join(real_data_dir, 'PRIME', 'capsid_protein_C.PRIME.json')
    with open(file_path, 'r') as f:
        return json.load(f)
