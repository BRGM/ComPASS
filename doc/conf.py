# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import sys
import os
import sphinx_rtd_theme

# This is for autodoc
sys.path.insert(0, os.path.abspath(".."))

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
    try:
        version_num = sys.argv[sys.argv.index("-t") + 1]
        if int(version_num[-1]) > 4:
            version = "5.x.x"
            scm_version = version
    except:
        pass
    return version, scm_version


version, release = get_version_info()

root_doc = "index"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.ifconfig",
    "recommonmark",  # for MarkDown
    "sphinx_revealjs",
    "sphinxcontrib.tikz",  # tikz pictures
    "scope",  #  to use meta and scope directive (skip file sometimes)
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "revealjs",
    "training/test",
    "training/meshes",
    "README.txt",
]

# to use custom config in .. ifconfig
def setup(app):
    app.add_config_value("versionlevel", "4", "env")


# -- Options for reveal.js output --------------------------------------------
# sphinx-build -M revealjs . .
revealjs_style_theme = "beige"
# contains static path (to store images or style files)
revealjs_static_path = ["_static"]

revealjs_script_conf = {
    "controls": True,
    "progress": True,
    "center": True,
    "transition": "slide",
}
revealjs_script_plugins = [
    {
        "name": "RevealNotes",
        "src": "revealjs4/plugin/notes/notes.js",
    },
    {
        "name": "RevealHighlight",
        "src": "revealjs4/plugin/highlight/highlight.js",
    },
]
revealjs_css_files = [
    "revealjs4/plugin/highlight/zenburn.css",
]

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
# sphinx-build -b latex . ./latex/; cd latex; make; cd -
latex_documents = [
    (
        "training/latex_index",
        "training_ComPASS_initiation.tex",
        "ComPASS Initiation",
        "L. Beaude, S. Lopez, F. Smai and various contributors",
        "manual",
    ),
]
latex_logo = "_static/compass_logo.png"
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    "papersize": "a4paper",
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    "figure_align": "!htb",
    "extraclassoptions": "openany",
}
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"
# contains static path (to store images or style files)
html_static_path = ["_static"]
html_css_files = ["custom.css"]  # todo to improve Solution appearance
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

rst_epilog = """
.. include:: <s5defs.txt>
.. raw:: html

    <style> .red {color:red} </style>
    <style> .lime {color:lime} </style>
"""
