# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "seismic-waveform-factory"
copyright = "2025, Thomas Ulrich"
author = "Thomas Ulrich"

release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"

# html_static_path = ['_static']

# If using autodoc for your Python code, set the main module here:
# autodoc_mock_imports = []

# If your code is in src/, add src to sys.path as well:
# sys.path.insert(0, os.path.abspath('../src'))
