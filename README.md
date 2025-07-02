# DMS Designer

A comprehensive toolkit for designing and analyzing deep mutational scanning (DMS) libraries.

## Overview

This repository contains tools for creating oligonucleotide libraries for protein mutagenesis studies, including:

- **Library Maker**: Generate comprehensive mutagenesis libraries with single and multiple amino acid mutants
- **Library Analysis**: Analyze and visualize DMS library composition and results

## Modules

### Library Maker
The `Library_maker` module provides functionality to create oligonucleotide libraries for deep mutational scanning experiments. It includes:

- Single and multiple amino acid mutant generation
- Synonymous nucleotide variant creation using codon usage tables
- Focused mutagenesis at specific positions
- Comprehensive library analysis and visualization

### Library Analysis  
The `library_analysis` module provides tools for analyzing DMS library results, including:

- Library composition analysis
- Variant distribution visualization
- Statistical analysis of mutation patterns
- Integration with MAVE-NN for deep learning analysis

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