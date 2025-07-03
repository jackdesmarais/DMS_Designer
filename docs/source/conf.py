# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import DMS_designer
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DMS Designer'
copyright = '2025, John J Desmarais'
author = 'John J Desmarais'
release = 'v0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'myst_parser',
    'nbsphinx',
]

# Source parsers and suffixes
source_suffix = ['.rst', '.md']

# Paths setup
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# MyST Parser settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
    'show-inheritance': True,
    'inherited-members': True
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Make sure the master doc is index
master_doc = 'index'

# Additional MyST Parser settings
myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Enable parsing of all files with .md extension
source_parsers = {
    '.md': 'myst_parser.sphinx_'
}

# Explicitly exclude setup.py from autodoc
autodoc_exclude_modules = ['setup']

# Add type hints
autodoc_typehints = 'both'
autodoc_typehints_description_target = 'documented'
autodoc_type_aliases = {}

# Show class and module docstrings
autoclass_content = 'both'

# Configure nbsphinx
nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
nbsphinx_kernel_name = 'python3'
nbsphinx_timeout = 600

# Add nbsphinx-specific settings
nbsphinx_prolog = """
{% set docname = env.doc2path(env.docname, base=None) %}

.. note::
   This page was generated from a Jupyter notebook.
   :download:`Download the notebook <{{ docname }}>`
"""