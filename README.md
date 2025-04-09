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
- **CFEL**: Contrast-FEL for comparing selection between clades
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
# Run analysis on HyPhy results
hyphy-results analyze --results-path /path/to/results --output-dir /path/to/output --site-mappings /path/to/mappings
```

### Input File Structure

The toolkit expects HyPhy results to be organized in the following structure:

```
results/
├── concat/
│   ├── BUSTED/
│   │   ├── gene1.BUSTED.json
│   │   ├── gene2.BUSTED.json
│   │   └── ...
│   ├── RELAX/
│   │   ├── gene1.RELAX.json
│   │   ├── gene2.RELAX.json
│   │   └── ...
│   └── contrastFEL/
│       ├── gene1.CFEL.json
│       ├── gene2.CFEL.json
│       └── ...
└── site_mappings/  # Optional, for coordinate mapping
    ├── gene1.csv
    ├── gene2.csv
    └── ...
```

#### Site Mapping Files (Optional)
CSV files mapping between different coordinate systems (e.g., between gene and protein coordinates). Format:
```csv
gene_position,protein_position,reference_codon
1,1,ATG
...
```

### Output Files
- **{gene}_summary.csv**: Gene-level summary statistics including:
  - RELAX: K parameter, p-values, and LRT statistics
  - BUSTED: Evidence of selection, omega distributions
  - Conservation metrics
  - Branch length information

- **{gene}_sites.csv**: Site-specific analysis including:
  - Beta values per clade
  - Substitution counts
  - Conservation status

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
