# hyphy-results-toolkit
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

## Usage

```bash
# Run analysis on HyPhy results from CAPHEINE workflow
hyphy-results -i /path/to/test_run/hyphy -o /path/to/output
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
├── PRIME/
│   ├── gene1.PRIME.json
│   └── ...
└── site_mappings/  # Optional, for coordinate mapping
    ├── gene1.csv
    ├── gene2.csv    
```

#### Site Mapping Files
TSV files mapping between site positions in different results (e.g., between FEL in clade1 and contrastFEL across all clades). Format:
```tsv
Consensus_Site	Consensus_capsid_protein_C_concat_1_Site	Consensus_capsid_protein_C_DENV1_1_Site
0	0	0
...
```

### Output Files
- **{gene}_summary.csv**: Gene-level summary statistics including:
  - RELAX: K parameter, p-values, and LRT statistics
  - BUSTED: Evidence of selection, omega distributions
  - Conservation metrics
  - Branch length information

- **{gene}_sites.csv**: Site-specific analysis including:
  - Beta values per comparison group
  - Substitution counts
  - Conservation status

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
pytest --cov=hyphy_results_toolkit tests/
```

The test suite includes:
- Unit tests for all HyPhy analysis methods
- Integration tests for gene processing
- Thread safety validation
- Mock data fixtures for reliable testing

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
