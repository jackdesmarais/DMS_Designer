"""
DMS Designer - Deep Mutational Scanning Library Design and Analysis Toolkit

A comprehensive toolkit for designing and analyzing deep mutational scanning (DMS) libraries.

This package provides tools for:
- Creating oligonucleotide libraries for protein mutagenesis studies
- Generating single and multiple amino acid mutants
- Creating synonymous nucleotide variants using codon usage tables
- Focused mutagenesis at specific positions
- Analyzing and visualizing library composition and results

Main Components:
- LibraryMaker: Main class for creating mutagenesis libraries
- LibraryAnalyzer: Tools for analyzing DMS library results
- Utility functions for codon/protein sequence conversion

Example Usage:
    >>> from DMS_designer import LibraryMaker
    >>> lm = LibraryMaker("ATGGCCGAA", single_mutants=True, doubles_to_make='equal')
    >>> 
    >>> from DMS_designer import summarize_file
    >>> summarize_file('my_library.csv', out_path='./analysis/')
"""

__version__ = "0.1.0"
__author__ = "John Desmarais"
__email__ = "desmara@cshl.edu"

# Import main classes and functions directly from the DMS_designer folder
try:
    from .library_maker import (
        LibraryMaker,
        aa_to_codon,
        protein_to_orf,
        orf_to_protein,
        summarize_file
    )
except ImportError as e:
    print(f"Warning: Could not import library_maker module: {e}")

try:
    from .library_analyzer import (
        Library,
        count_over_thresh,
        make_density,
        line_plotter,
        correspondance_plotter,
        density_scatter,
        gini,
        multi_gini,
        depth_plot,
        make_ROC,
    )
except ImportError as e:
    print(f"Warning: Could not import library_analyzer module: {e}")

# Define what gets imported with "from DMS_designer import *"
__all__ = [
    'LibraryMaker',
    'Library', 
    'summarize_file'
]

# Package metadata
__package_info__ = {
    "name": "DMS_designer",
    "version": __version__,
    "description": "Deep Mutational Scanning Library Design and Analysis Toolkit",
    "author": __author__,
    "email": __email__,
    "modules": [
        "library_maker",
        "library_analyzer"
    ]
}
