# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "KN1D"
copyright = "2026, Jamie Dunsmore"
author = "Jamie Dunsmore"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.coverage",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Numpy-doc config
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_attr_annotations = True
napoleon_preprocess_types = True

# Enable numbered references
numfig = True

autodoc_type_aliases = {
    "NDArray": "numpy.typing.NDArray",
    "ArrayLike": "numpy.typing.ArrayLike",
}

napoleon_type_aliases = {
    "ndarray": "numpy.typing.NDArray",
}

autodoc_default_options = {"ignore-module-all": True}
autodoc_typehints = "description"
# autodoc_class_signature = "mixed"

# The default role for text marked up `like this`
default_role = "any"

# Tell sphinx what the primary language being documented is.
primary_domain = "py"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "python"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/W-M-plasma-group/KN1DPy",
}

pygments_style = "sphinx"

# -- Options for intersphinx extension ---------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}

nitpick_ignore = {
    ("py:class", "numpy.typing.ArrayLike"),
    ("py:class", "numpy.typing.NDArray"),
    ("py:class", "numpy.typing.Scalar"),
}
nitpick_ignore_regex = {
    ("py:class", r".*_ScalarT"),
}
