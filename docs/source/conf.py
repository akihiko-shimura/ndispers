# Configuration file for the Sphinx documentation builder.

import ndispers

# -- Project information

project = 'ndispers'
copyright = '2021, Akihiko Shimura'
author = 'Akihiko Shimura'

version = ndispers.__version__
release = version

# -- General configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# The class docstrings use informal section titles; napoleon must know them
# or Sphinx rejects the rST section underlines inside docstrings.
napoleon_custom_sections = [
    'Sellmeier equation',
    'Thermo-optic coefficient',
    'Validity range',
    'Ref',
    'Usage',
    'Note',
    'Example',
    'input',
    'return',
]

autodoc_member_order = 'bysource'

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
