# DMS Designer

A toolkit for designing and analyzing deep mutational scanning (DMS) libraries. Associated with the paper Alpsoy, Ipsaro, Skopelitis, et al. Structural Basis of DNA-Dependent Coactivator Recruitment by the Tuft Cell Master Regulator POU2F3.

## Overview

This repository contains tools for creating libraries of oligonucleotides for protein mutagenesis studies as well as tools for analyzing the results of such studies, including:

- **Library Maker**: Generate comprehensive mutagenesis libraries with single and multiple amino acid mutants
- **Library Analysis**: Analyze and visualize mutagenesis library results

## Modules

### Library Maker
The `Library_maker` module provides functionality to create oligonucleotide libraries for deep mutational scanning experiments. It includes:

- Single and multiple amino acid mutant generation
- Synonymous nucleotide variant creation using codon usage tables
- Focused mutagenesis at specific positions
- Comprehensive library analysis and visualization

### Library Analysis  
The `library_analysis` module provides tools for analyzing DMS library results, including:

- Creation of quality control plots
   - Measuring library composition and skew
   - Comparing replicates
   - Assessing selection at time point 0
- library scale correction
   - Replicate scale correction for merging
   - Cross library scale correction for merging
- Variant fitness calculation and plotting
   - Create variant effect scores from library counts
   - Create single mutant variant effect heatmaps

## Installation

1. Clone this repository
2. Install dependencies using the provided environment files:
   ```bash
   conda env create -f env.yml
   ```

## Usage

### Creating a Library
```python
from Library_maker.library_maker import LibraryMaker

# Create a simple library with single mutants
lm = LibraryMaker("ATGGCCGAA", single_mutants=True, doubles_to_make=0)
```

### Analyzing a Library
```python
from library_analysis.library_analyzer import LibraryAnalyzer

# Analyze library composition
analyzer = LibraryAnalyzer("library_file.csv")
analyzer.generate_summary_plots()
```

## Documentation

For detailed documentation, see the docstrings in the individual modules or run:
```python
help(LibraryMaker)
```

## License

[Add your license information here] 