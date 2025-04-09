"""
Pytest configuration and fixtures.
"""

import os
import json
import tempfile
from typing import Dict, Any
import pytest

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield tmpdirname

@pytest.fixture
def mock_results_dir(temp_dir):
    """Create a mock results directory with HyPhy result files."""
    # Create directory structure
    for method in ['BUSTED', 'RELAX', 'CFEL', 'FEL', 'MEME', 'PRIME']:
        os.makedirs(os.path.join(temp_dir, 'concat', method), exist_ok=True)
    
    yield temp_dir

@pytest.fixture
def mock_busted_results() -> Dict[str, Any]:
    """Create mock BUSTED results."""
    return {
        'test results': {
            'p-value': 0.01,
            'LRT': 10.5
        },
        'fits': {
            'Unconstrained model': {
                'Rate Distributions': {
                    'Test': {
                        '2': {'omega': 3.5, 'weight': 0.1}
                    }
                }
            }
        },
        'branch attributes': {
            '0': {
                'branch1': {'length': 0.1},
                'branch2': {'length': 0.2}
            }
        }
    }

@pytest.fixture
def mock_relax_results() -> Dict[str, Any]:
    """Create mock RELAX results."""
    return {
        'test results': {
            'p-value': 0.05,
            'relaxation or intensification parameter': 1.5,
            'LRT': 8.2
        }
    }

@pytest.fixture
def mock_cfel_results() -> Dict[str, Any]:
    """Create mock CFEL results."""
    return {
        'tested': {
            '0': {
                'branch1': 'clade1',
                'branch2': 'clade2'
            }
        },
        'MLE': {
            'headers': [
                ('site', 'Site'),
                ('beta (clade1)', 'Beta for clade1'),
                ('subs (clade1)', 'Substitutions for clade1'),
                ('beta (clade2)', 'Beta for clade2'),
                ('subs (clade2)', 'Substitutions for clade2')
            ],
            'content': {
                '0': [
                    [1, 0.5, 2, 1.0, 0],
                    [2, 1.0, 0, 1.0, 0]
                ]
            }
        }
    }

@pytest.fixture
def write_mock_results(mock_results_dir, mock_busted_results, mock_relax_results, mock_cfel_results):
    """Write mock results to files."""
    def _write_results(gene: str):
        # Write BUSTED results
        with open(os.path.join(mock_results_dir, 'concat', 'BUSTED', f'{gene}.BUSTED.json'), 'w') as f:
            json.dump(mock_busted_results, f)
        
        # Write RELAX results
        with open(os.path.join(mock_results_dir, 'concat', 'RELAX', f'{gene}.RELAX.json'), 'w') as f:
            json.dump(mock_relax_results, f)
        
        # Write CFEL results
        with open(os.path.join(mock_results_dir, 'concat', 'contrastFEL', f'{gene}.CFEL.json'), 'w') as f:
            json.dump(mock_cfel_results, f)
    
    return _write_results
