# DRHIP: Data Reduction for HyPhy with Inference Processing
A Python toolkit for analyzing HyPhy results. Process and summarize evolutionary selection analysis results across multiple methods and genes.

## Installation

```bash
pip install -e .
```

## Features

### Analysis Methods Support
- **RELAX**: Test for relaxation of selection
- **BUSTED**: Branch-site unrestricted statistical test for episodic diversification
- **CFEL**: Contrast-FEL for comparing selection between groups
- **FEL**: Fixed effects likelihood test for selection
- **MEME**: Mixed effects model of evolution
- **PRIME**: Property Informed Models of Evolution

### Key Features
- Process HyPhy results for multiple genes in parallel
- Generate comprehensive summary statistics
- Site-specific analysis with conservation tracking
- Thread-safe output handling
- **Combined output files** across all genes for easy analysis

## Usage

```bash
# Run analysis on HyPhy results from CAPHEINE workflow
drhip -i /path/to/results/hyphy -o /path/to/output
```

### Input File Structure

The toolkit is designed to work with HyPhy results from the CAPHEINE workflow, which are organized in the following structure:

```
hyphy/
├── BUSTED/
│   ├── gene1.BUSTED.json
│   └── ...
├── FEL/
│   ├── gene1.FEL.json
│   └── ...
├── MEME/
│   ├── gene1.MEME.json
│   └── ...
└── PRIME/
    ├── gene1.PRIME.json
    └── ...
```

### Output Files

The toolkit produces combined files that aggregate data across all genes:

- **combined_summary.csv**: Gene-level summary statistics including:
  - RELAX: K parameter, p-values, and LRT statistics
  - BUSTED: Evidence of selection, omega distributions
  - Conservation metrics
  - Branch length information

- **combined_sites.csv**: Site-specific analysis including:
  - Beta values per comparison group
  - Substitution counts
  - Conservation status

- **combined_comparison_summary.csv**: Comparison group-specific summary statistics (when comparison groups are present):
  - Group sizes (N)
  - Total branch lengths (T)
  - Group-specific dN/dS ratios
  - Conservation metrics per group

- **combined_comparison_site.csv**: Comparison group-specific site analysis (when comparison groups are present):
  - Site-specific metrics for each comparison group
  - Conservation markers
  - Group-specific selection indicators

The combined files automatically include the superset of all fields found across all genes, with missing values marked as 'NA'. This makes it easy to analyze patterns across the entire dataset.

### Data Representation

The toolkit uses consistent markers to represent different types of data:

- **Significant Results**: Actual values (p-values, property names, etc.)
- **Non-significant Results**: "-" marker for most fields
- **Missing or Malformed Data**: "NA" marker

If data for a particular method is missing for a gene, or if the data is malformed or cannot be processed, the toolkit will output "NA" in the corresponding fields. This allows for easy identification of missing data versus non-significant results.

## Extending the Toolkit

### Adding New Fields for Specific Methods

To add new fields to the output for an existing method:

1. **Modify the method's class**: Each method has its own class in `drhip/methods/` (e.g., `busted.py`, `relax.py`).
   - Add the new field to the `process_results()` method for summary data
   - Add the new field to the `process_site_data()` method for site-specific data

2. **Update validation fields** (optional): If the field should be validated as a required field:
   - Add it to `SUMMARY_FIELDNAMES` or `SITES_FIELDNAMES` in `config.py`

Example for adding a new summary field to BUSTED:

```python
def process_results(self, results):
    summary_data = super().process_results(results)
    
    # Add your new field
    summary_data['new_field_name'] = self._extract_new_field(results)
    
    return summary_data
    
def _extract_new_field(self, results):
    # Logic to extract the new field from results
    # Return 'NA' if data is missing or invalid
    try:
        # Extract data from results
        return extracted_value
    except:
        return 'NA'
```

### Field Validation

The toolkit validates that all expected fields are present in the output:

1. **Required fields** are defined in `config.py`:
   - `SUMMARY_FIELDNAMES`: Required fields for summary files
   - `SITES_FIELDNAMES`: Required fields for site files

2. **Field validation** occurs in `process_gene.py`:
   - Missing fields are reported as warnings
   - This ensures consistent output structure across all genes

### Adding a New Analysis Method

To add support for a new HyPhy analysis method:

1. **Create a new method class** in `drhip/methods/`:
   - Inherit from `HyPhyMethod` in `base.py`
   - Implement the required methods: `process_results()` and optionally `process_site_data()`

2. **Register the method** in `registry.py`:
   - Import your new method class
   - Add it to the `__init__` method of `HyPhyMethodRegistry`

3. **Update the method paths** in `config.py`:
   - Add an entry to `METHOD_PATHS` if your method uses a different directory name

Example for a new method:

```python
# new_method.py
from .base import HyPhyMethod

class NewMethod(HyPhyMethod):
    def __init__(self):
        super().__init__(name='NEW_METHOD', file_suffix='NEW_METHOD.json')
        
    def process_results(self, results):
        # Process and return summary data
        return {...}
        
    def process_site_data(self, results):
        # Process and return site-specific data
        return {...}
```

Then in `registry.py`:
```python
from .new_method import NewMethod
# ...
self.register(NewMethod())
```

### Comparison Groups Detection

Comparison groups represent different sets of branches in the phylogenetic tree (e.g., 'foreground' vs 'background', or 'test' vs 'reference') that are being compared for evolutionary selection differences:

1. **What are comparison groups**: 
   - Groups of branches in the phylogenetic tree labeled for comparative analysis
   - Used by methods like CFEL and RELAX to compare selection pressures between different lineages
   - Typically labeled as 'test'/'reference', 'foreground'/'background', or custom labels

2. **Detection logic** is in `utils/result_helpers.py`:
   - `detect_comparison_groups()` attempts to extract group labels from CFEL and RELAX results files
   - The toolkit examines these results to find branch labels assigned to different groups
   - Falls back to default groups ('test'/'reference') if no explicit labels are detected

3. **Adding detection for a new method**:
   - Add a detection function in `detect_comparison_groups()`
   - Add your method to the `detection_functions` dictionary
   - Add your method to the `methods_to_check` list

4. **Consistency validation**:
   - The toolkit ensures that groups are consistent across methods
   - Raises an error if inconsistent groups are detected

5. **Method-specific requirements**:
   - Some methods (like RELAX) require comparison groups to be set before returning fields
   - The toolkit gracefully handles these requirements with proper error handling
   - If a method can't provide fields without comparison groups, those fields are skipped

6. **CSV combining logic**:
   - The toolkit combines CSV files with similar suffixes across all genes
   - Handles various output types: summary, sites, comparison_summary, comparison_site, and comparison
   - Ensures proper field ordering for each output type

## Development and Testing

### Setting up Development Environment

Install the package with development dependencies:

```bash
pip install -e ".[dev]"
```

### Running Tests

To run the test suite:

```bash
pytest tests/
```

To run tests with coverage reporting:

```bash
pytest --cov=drhip tests/
```

The test suite includes:
- Unit tests for all HyPhy analysis methods
- Integration tests for gene processing
- Thread safety validation
- Minimal actual HyPhy output files for reliable testing

## Requirements

- Python 3.7+
- numpy
- pandas
- scipy

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authors

- Sergei L Kosakovsky Pond (spond@temple.edu)
- Hannah Verdonk (hannah.verdonk@temple.edu)
- Danielle Callan (dcallan@temple.edu)
