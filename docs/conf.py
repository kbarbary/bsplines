import sys
import os
import sphinx_rtd_theme
import matplotlib

# generate api directory if it doesn't already exist
if not os.path.exists('api'):
    os.mkdir('api')

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'numpydoc',
    'matplotlib.sphinxext.plot_directive']

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

master_doc = 'index'
numpydoc_show_class_members = False
autosummary_generate = ["index.rst"]
autoclass_content = "class"
autodoc_default_flags = ["members", "no-special-members"]

project = 'bsplines'
author = 'Kyle Barbary and contributors'
copyright = '2016, ' + author
version = "0.1"
release = "0.1.0"
