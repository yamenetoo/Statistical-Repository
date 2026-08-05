import os
import sys
sys.path.insert(0, os.path.abspath('../src'))

project = 'unit_transform_dist'
copyright = '2026, Mohamad Yamen AL Mohamad, Ali Akbar Heydari, Hossein Jabbari Khamnei'
author = 'Mohamad Yamen AL Mohamad, Ali Akbar Heydari, Hossein Jabbari Khamnei'
release = '0.3.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
