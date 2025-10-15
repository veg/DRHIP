"""
Tests for input JSON validation and orchestrator skip behavior.
"""

import os
import json
import shutil
import tempfile
import csv
import pytest

from drhip.methods import (
    BustedMethod,
    FelMethod,
    CfelMethod,
)
from drhip.parsers.process_gene import process_gene


def test_busted_validate_input_json_reports_missing():
    method = BustedMethod()
    # Empty results should miss both required paths
    results = {}
    missing = set(method.validate_input_json(results))
    assert 'test results.p-value' in missing
    assert 'fits.Unconstrained model.Rate Distributions.Test' in missing
    assert 'input.trees.0' in missing
    assert 'substitutions.0' in missing

def test_busted_validate_input_json_works_correctly():
    method = BustedMethod()
    # Present MLE headers/content but no required header names
    results = {
        'test results': {
            'p-value': 0.0001,
        },
        'fits': {
            'Unconstrained model': {
                'Rate Distributions': {
                    'Test': 0.0001
                }
            }
        },
         "input":{
            "trees":{
                "0":"(PP563828_1_2023_09_04,PP563831_1_2023_08_29,PP563880_1_2023_10_05,((PP564823_1_2023_10_06,PP563845_1_2023_10_22)Node4,(PP563884_1_2023_09_24,PP563838_1_2023_09_30,PP563839_1_2023_09_29)Node7)Node3,PP563832_1_2023_08_28)"
                }
            },
        "substitutions":{
            "0":{
                "0":{
                "root":"ATG"
                },
                "1":{
                "root":"AAC"
                },
            }
        }
    }
    missing = set(method.validate_input_json(results))
    assert len(missing) == 0


def test_fel_validate_input_json_requires_headers_and_columns():
    method = FelMethod()
    # Present MLE headers/content but no required header names
    results = {
        'MLE': {
            'headers': [["p-value", "Asymptotic p-value for evidence of selection, i.e. beta &neq; alpha"],],
            'content': {'0': []}
        }
    }
    missing = set(method.validate_input_json(results))
    # Should report missing header names
    assert 'MLE.headers:alpha' in missing
    assert 'MLE.headers:beta' in missing
    assert 'MLE.headers:p-value' not in missing


def test_cfel_validate_input_json_requires_mle_and_tested():
    method = CfelMethod()
    # Missing substitutions data (hyphy <=v2.5.83)
    results = {
        'MLE': {
            'headers': [
                ["alpha", "Synonymous substitution rate at a site"],
                ["beta (background)", "Non-synonymous substitution rate at a site for background branches"],
                ["beta (Reference)", "Non-synonymous substitution rate at a site for Reference branches"],
                ["beta (Foreground)", "Non-synonymous substitution rate at a site for Foreground branches"],
                ["subs (Reference)", "Substitutions mapped to Reference branches"],
                ["subs (Foreground)", "Substitutions mapped to Foreground branches"],
                ["P-value (overall)", "P-value for the test that non-synonymous rates differ between any of the selected groups: Reference,Foreground"],
                ["Q-value (overall)", "Q-value for the test that non-synonymous rates differ between any of the selected groups: Reference,Foreground"],
                ["Permutation p-value", "Label permutation test for significant sites"],
                ["Total branch length", "The total length of branches contributing to inference at this site, and used to scale beta-alpha"] 
            ],
            'content': {'0': [0, 0, 0, 0, 1, 1, 1, 1, 1, 0]}
        },
        "tested":{
            "0":{
                "AB178040_1_2002":"background",
                "AB195673_1_2003":"background",
            }
        },
        "branch attributes":{
            "0":{
                "AB178040_1_2002":{
                "Global MG94xREV":0.003354111713152422,
                "Nucleotide GTR":0.003345420161232569,
                "original name":"AB178040_1_2002"
                },
            },
        },
        "fits":{
            "Global MG94xREV":{
                "Rate Distributions":{
                    "Test":0.0001
                }
            }
        },
        "input":{
            "trees":{
                "0":"(PP564823_1_2023_10_06,((((((AY732475_1_1994,AY732480_1_1994,(((PP563832_1_2023_08_28,PP563826_1_2023_08_21,PP563831_1_2023_08_29)Node12,(PP563838_1_2023_09_30,PP563839_1_2023_09_29)Node17)Node11,PP563845_1_2023_10_22)Node10,(AY732482_1_2001,AB178040_1_2002)Node21)Node7,AY726553_1_2002)Node6,AY708047_1_2001)Node5,((AY732478_1_1991,AY726552_1_2002,AY732477_1_1991)Node29,AF298808_1_1998)Node28)Node4,(AB204803_1_2004,AB195673_1_2003)Node35)Node3,AF298807_1_1998)Node2,AY732474_1_1980)"
                }
            },
    }
    missing = set(method.validate_input_json(results))
    assert 'tested.0' not in missing
    assert 'branch attributes.0.Global MG94xREV' not in missing
    assert 'fits.Global MG94xREV.Rate Distributions' not in missing
    assert 'substitutions.0' in missing
    assert len(missing) == 1  # WRONG once all things are added as required, but it should just be subs missing


def test_process_gene_warns_invalid_busted(results_dir, capsys):
    
    gene_name = 'capsid_protein_C'

    # Use a temp copy of the real results, but replace BUSTED JSON for the gene 
    # with a fake json that's missing the "substitutions" fields
    results = {
        "tested":{
            "0":{
                "AB178040_1_2002":"background",
                "AB195673_1_2003":"background",
            }
        },
        "branch attributes":{
            "0":{
                "AB178040_1_2002":{
                "Global MG94xREV":0.003354111713152422,
                "Nucleotide GTR":0.003345420161232569,
                "original name":"AB178040_1_2002"
                },
            },
        },
        "fits":{
            'Unconstrained model': {
                'Rate Distributions': {
                    "Test":{
                        "0":{
                        "omega":226.0589415365084,
                        "proportion":0
                        },
                        "1":{
                        "omega":0.04198616898833647,
                        "proportion":0.9147510556650086
                        },
                        "2":{
                        "omega":0.1430120656070995,
                        "proportion":0.08477641331245425
                        },
                        "3":{
                        "omega":1,
                        "proportion":0.0004725310225371633
                        }
                    }
                }
            }
        },
        "test results":{
            "LRT":0,
            "p-value":0.5
        },
        "input":{
            "trees":{
                "0":"(PP564823_1_2023_10_06,((((((AY732475_1_1994,AY732480_1_1994,(((PP563832_1_2023_08_28,PP563826_1_2023_08_21,PP563831_1_2023_08_29)Node12,(PP563838_1_2023_09_30,PP563839_1_2023_09_29)Node17)Node11,PP563845_1_2023_10_22)Node10,(AY732482_1_2001,AB178040_1_2002)Node21)Node7,AY726553_1_2002)Node6,AY708047_1_2001)Node5,((AY732478_1_1991,AY726552_1_2002,AY732477_1_1991)Node29,AF298808_1_1998)Node28)Node4,(AB204803_1_2004,AB195673_1_2003)Node35)Node3,AF298807_1_1998)Node2,AY732474_1_1980)"
                }
            },
    }

    with tempfile.TemporaryDirectory() as tmp_results:
        # Copy entire results tree
        for entry in os.listdir(results_dir):
            src_path = os.path.join(results_dir, entry)
            dst_path = os.path.join(tmp_results, entry)
            if os.path.isdir(src_path):
                shutil.copytree(src_path, dst_path)

        # Overwrite the gene's BUSTED file with invalid JSON
        busted_dir = os.path.join(tmp_results, 'BUSTED')
        os.makedirs(busted_dir, exist_ok=True)
        busted_file = os.path.join(busted_dir, f'{gene_name}.BUSTED.json')
        with open(busted_file, 'w') as f:
            json.dump(results, f)

        with tempfile.TemporaryDirectory() as output_dir:
            process_gene(gene_name, tmp_results, output_dir)

            # Capture logs and ensure BUSTED was skipped due to validation
            captured = capsys.readouterr()
            assert "WARNING: BUSTED for capsid_protein_C is missing required fields: ['substitutions.0']." in captured.out
            assert gene_name in captured.out

            # Validate summary output exists and contains BUSTED and non-BUSTED fields
            summary_file = os.path.join(output_dir, f'{gene_name}_summary.csv')
            assert os.path.exists(summary_file)
            with open(summary_file, 'r') as f:
                reader = csv.DictReader(f)
                print(row)
                # Ensure that present BUSTED fields are present in the output
                assert 'BUSTED_pval' in row
                assert 'BUSTED_omega3' in row
                # FEL/MEME fields should still be present
                assert 'negative_sites' in row
                assert 'N' in row

            # Validate sites output exists and contains non-BUSTED fields, but not BUSTED substitution data
            sites_file = os.path.join(output_dir, f'{gene_name}_sites.csv')
            assert os.path.exists(sites_file)
            with open(sites_file, 'r') as f:
                reader = csv.DictReader(f)
                row = next(reader)
                # Ensure that absent BUSTED fields are absent in the output
                assert 'composition' not in row
                assert 'substitutions' not in row
                # FEL/PRIME fields should still be present
                assert 'fel_selection' in row
                assert 'prime_marker' in row
