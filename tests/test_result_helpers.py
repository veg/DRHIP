"""
Tests for result helper functions.
"""

import pytest
from unittest.mock import patch, MagicMock

from hyphy_results_toolkit.utils import result_helpers as rh


def test_merge_field_value_new_field():
    """Test merging a new field value."""
    data_dict = {}
    providers = {}
    
    rh.merge_field_value(data_dict, "test_field", "test_value", "test_method", providers)
    
    assert data_dict["test_field"] == "test_value"
    assert providers["test_field"] == ["test_method"]


def test_merge_field_value_same_value():
    """Test merging a field with the same value."""
    data_dict = {"test_field": "test_value"}
    providers = {"test_field": ["method1"]}
    
    rh.merge_field_value(data_dict, "test_field", "test_value", "method2", providers)
    
    assert data_dict["test_field"] == "test_value"
    assert providers["test_field"] == ["method1", "method2"]


def test_merge_field_value_na_existing():
    """Test merging when existing value is 'NA'."""
    data_dict = {"test_field": "NA"}
    providers = {"test_field": ["method1"]}
    
    rh.merge_field_value(data_dict, "test_field", "test_value", "method2", providers)
    
    assert data_dict["test_field"] == "test_value"
    assert providers["test_field"] == ["method1", "method2"]


def test_merge_field_value_na_new():
    """Test merging when new value is 'NA'."""
    data_dict = {"test_field": "test_value"}
    providers = {"test_field": ["method1"]}
    
    rh.merge_field_value(data_dict, "test_field", "NA", "method2", providers)
    
    assert data_dict["test_field"] == "test_value"
    assert providers["test_field"] == ["method1", "method2"]


def test_merge_field_value_conflict():
    """Test merging when values conflict."""
    data_dict = {"test_field": "value1"}
    providers = {"test_field": ["method1"]}
    
    with patch('builtins.print') as mock_print:
        rh.merge_field_value(data_dict, "test_field", "value2", "method2", providers)
        
        # Verify warning was printed
        assert mock_print.call_count >= 1
        
        # Verify new value was used
        assert data_dict["test_field"] == "value2"
        assert providers["test_field"] == ["method1", "method2"]


def test_merge_field_value_with_context():
    """Test merging with context information."""
    data_dict = {"test_field": "value1"}
    providers = {"test_field": ["method1"]}
    context = {"gene": "gene1", "site": 42}
    
    with patch('builtins.print') as mock_print:
        rh.merge_field_value(data_dict, "test_field", "value2", "method2", providers, context)
        
        # Verify warning was printed with context
        assert mock_print.call_count >= 1
        # Check that context was included in at least one message
        context_included = False
        for call in mock_print.call_args_list:
            args = call[0][0]
            if isinstance(args, str) and "gene1" in args and "42" in args:
                context_included = True
                break
        assert context_included


def test_merge_method_data():
    """Test merging entire method data dictionaries."""
    target_dict = {"field1": "value1"}
    method_data = {"field2": "value2", "field3": "value3"}
    providers = {"field1": ["method1"]}
    
    rh.merge_method_data(target_dict, method_data, "method2", providers)
    
    assert target_dict == {
        "field1": "value1",
        "field2": "value2",
        "field3": "value3"
    }
    assert providers == {
        "field1": ["method1"],
        "field2": ["method2"],
        "field3": ["method2"]
    }


def test_ensure_ordered_fields():
    """Test ensuring priority fields come first."""
    all_fields = {"field1", "field2", "field3", "field4", "field5"}
    priority_fields = ["field3", "field1", "field6"]  # field6 not in all_fields
    
    ordered = rh.ensure_ordered_fields(all_fields, priority_fields)
    
    # Check that priority fields come first in the specified order
    assert ordered[:2] == ["field3", "field1"]
    # Check that all fields are included
    assert set(ordered) == all_fields
    # Check that non-priority fields are sorted
    assert ordered[2:] == sorted(["field2", "field4", "field5"])


def test_collect_method_fields():
    """Test collecting fields from methods with results."""
    # Create mock methods with proper code attribute handling
    def create_mock_method(name, fields):
        method = MagicMock()
        method.name = name
        
        # Create a real function to mock get_summary_fields
        def get_fields(results=None):
            return fields
        
        # Assign the real function to the mock
        method.get_summary_fields = get_fields
        return method
    
    method1 = create_mock_method("Method1", ["field1", "field2"])
    method2 = create_mock_method("Method2", ["field2", "field3"])
    method3 = create_mock_method("Method3", ["field4"])
    
    # Set up method results
    method_results = {
        "Method1": {"result": "data1"},
        "Method3": {"result": "data3"}
        # Method2 has no results
    }
    
    # Test collecting fields
    fields = rh.collect_method_fields(
        [method1, method2, method3],
        method_results,
        "get_summary_fields"
    )
    
    # Should only include fields from methods with results
    assert fields == {"field1", "field2", "field4"}


def test_validate_fields_all_present():
    """Test field validation when all expected fields are present."""
    expected_fields = {"field1", "field2"}
    output_fields = {"field1", "field2", "field3"}
    
    with patch('builtins.print') as mock_print:
        rh.validate_fields(expected_fields, output_fields)
        
        # No warning should be printed
        mock_print.assert_not_called()


def test_validate_fields_missing():
    """Test field validation when some expected fields are missing."""
    expected_fields = {"field1", "field2", "field3"}
    output_fields = {"field1", "field3"}
    
    with patch('builtins.print') as mock_print:
        rh.validate_fields(expected_fields, output_fields, "summary", "gene1")
        
        # Warning should be printed
        mock_print.assert_called()
        # Check that missing field and context are mentioned
        args = mock_print.call_args[0][0]
        assert "field2" in args
        assert "summary" in args
        assert "gene1" in args


def test_detect_comparison_groups():
    """Test detecting comparison groups from method results."""
    # Mock method results with comparison groups
    method_results = {
        "RELAX": {
            "test results": {
                "relaxation or intensification parameter": {
                    "foreground": 0.5,
                    "background": 1.0,
                    "overall": 0.75
                }
            }
        },
        "CFEL": {
            "tested": {
                "0": {
                    "branch1": "foreground",
                    "branch2": "background"
                }
            }
        }
    }
    
    # Test detection
    groups, groups_by_method = rh.detect_comparison_groups(method_results)
    
    # Verify results
    assert set(groups) == {"foreground", "background"}
    
    # Check that the groups are detected from both methods
    assert "RELAX" in groups_by_method
    assert "CFEL" in groups_by_method
    
    # Check that both methods detected the same groups (order doesn't matter)
    assert set(groups_by_method["RELAX"]) == {"foreground", "background"}
    assert set(groups_by_method["CFEL"]) == {"foreground", "background"}


def test_detect_comparison_groups_default():
    """Test detecting comparison groups with default values."""
    # Empty method results
    method_results = {}
    default_groups = ["group1", "group2"]
    
    # Test detection with defaults
    groups, groups_by_method = rh.detect_comparison_groups(method_results, default_groups)
    
    # Verify results
    assert groups == default_groups
    assert groups_by_method == {}


def test_detect_comparison_groups_inconsistent():
    """Test detecting inconsistent comparison groups."""
    # Mock method results with inconsistent groups
    method_results = {
        "RELAX": {
            "test results": {
                "relaxation or intensification parameter": {
                    "foreground": 0.5,
                    "background": 1.0,
                    "overall": 0.75
                }
            }
        },
        "CFEL": {
            "tested": {
                "0": {
                    "branch1": "group1",
                    "branch2": "group2"
                }
            }
        }
    }
    
    # Test detection with inconsistent groups
    with pytest.raises(ValueError):
        rh.detect_comparison_groups(method_results)
