import sys
import os
import sphinx_rtd_theme
import matplotlib

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'matplotlib.sphinxext.plot_directive']

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

master_doc = 'index'
autosummary_generate = True
autoclass_content = "class"
autodoc_default_flags = ["members", "no-special-members"]

project = 'bsplines'
author = 'Kyle Barbary and contributors'
copyright = '2016, ' + author
version = "0.1"
release = "0.1.0"
