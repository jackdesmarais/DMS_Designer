#!/usr/bin/env python3
"""
Setup script for DMS_designer package.

Deep Mutational Scanning Library Design and Analysis Toolkit
"""

from setuptools import setup, find_packages
import os
import re

# Read the README file for long description
def read_readme():
    """Read README.md file for long description."""
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

# Read version from __init__.py
def get_version():
    """Extract version from __init__.py file."""
    init_path = os.path.join(os.path.dirname(__file__), 'DMS_designer', '__init__.py')
    with open(init_path, 'r', encoding='utf-8') as f:
        content = f.read()
        version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
        if version_match:
            return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

# Package configuration
setup(
    name="DMS_designer",
    version=get_version(),
    author="John Desmarais",
    author_email="desmara@cshl.edu",
    description="Deep Mutational Scanning Library Design and Analysis Toolkit",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/jackdesmarais/DMS_Designer",
    project_urls={
        "Bug Reports": "https://github.com/jackdesmarais/DMS_Designer/issues",
        "Source": "https://github.com/jackdesmarais/DMS_Designer",
        "Documentation": "https://github.com/jackdesmarais/DMS_Designer#readme",
    },
    packages=find_packages(),
    package_data={
        'DMS_designer': ['codon_usage.xlsx'],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "logomaker>=0.8.0",
        "openpyxl>=3.0.0",
        "tqdm>=4.60.0",
        "biopython>=1.79.0",
        "mavenn",
        "statsmodels",
    ],
    extras_require={
        "docs": [
            "sphinx==5.0.2",
            "sphinx_rtd_theme",
            "myst-parser",
            "nbsphinx",  # For rendering Jupyter notebooks
            "pandoc", # For rendering Jupyter notebooks
            "ipython",   # Required by nbsphinx
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "dms-maker=DMS_designer.library_maker:main",
        ],
    },
    keywords=[
        "deep mutational scanning",
        "protein mutagenesis",
        "oligonucleotide library",
        "bioinformatics",
        "protein engineering",
        "DMS",
        "mutagenesis",
        "library design",
    ],
    zip_safe=False,
    platforms=["any"],
    license="MIT",
) 