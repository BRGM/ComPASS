# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import sys
import os
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = "ComPASS"
copyright = "2013-2019, various contributors"
author = "various contributors"


def get_version_info():
    import re
    from setuptools_scm import get_version

    scm_version = get_version(root="..", relative_to=__file__)
    matches = re.match("(\d+\.\d+.\d+)(.*)", scm_version)
    version, tag = matches.groups()
    if len(tag) > 0:
        assert tag.startswith(".dev")
        matches = re.match("(\d+\.\d+)(.*)", version)
        major = matches.group(1)
        version = f"{major}.x"
    return version, scm_version


version, release = get_version_info()

master_doc = "index"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "recommonmark",  # for MarkDown
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"
# contains static path (to store images or style files)
html_static_path = ["_static"]

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "_static/compass_logo.png"

# Add a timestamp using datetime format
html_last_updated_fmt = "%b %d, %Y, %X"

# The suffix(es) of source filenames.
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}
