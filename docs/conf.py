"""Sphinx configuration: ``uv run sphinx-build -b html docs docs/_build/html``."""

import tomllib
from pathlib import Path

from ferrospot import __version__

_meta = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())

project = "ferrospot"
author = ", ".join(a["name"] for a in _meta["project"]["authors"])
copyright = f"2026, {author}"
version = release = __version__

extensions = [
    "myst_parser",
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "**/README.md"]

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "amsmath", "tasklist"]
myst_heading_anchors = 3

# autoapi documents src/ferrospot from the docstrings (NumPy style).
autoapi_dirs = ["../src/ferrospot"]
autoapi_options = ["members", "undoc-members", "show-module-summary"]
autoapi_ignore = ["*dashboard*"]
autoapi_add_toctree_entry = True
autodoc_typehints = "none"
always_document_param_types = True
typehints_fully_qualified = False
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True
suppress_warnings = ["ref.python", "ref.class", "ref.mod", "ref.func"]

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
}

html_theme = "shibuya"
html_title = "ferrospot"
html_static_path = ["_static"]
html_theme_options = {
    "accent_color": "orange",
    "github_url": "https://github.com/arunoruto/ferrospot",
}
# "Edit this page" and repository links (docs in /docs/ on main).
html_context = {
    "source_type": "github",
    "source_user": "arunoruto",
    "source_repo": "ferrospot",
}


# Pages that stay local (git-ignored: results on copyrighted figures) get a
# toctree on the index only where they exist, so CI builds without them.
LOCAL_PAGES = [
    p
    for p in ("regression", "calibration_results")
    if (Path(__file__).parent / f"{p}.md").exists()
]


def _local_toctree(app, docname, source):
    if docname == "index" and LOCAL_PAGES:
        entries = "\n".join(LOCAL_PAGES)
        source[0] += (
            f"\n```{{toctree}}\n:maxdepth: 1\n:caption: Local results\n\n{entries}\n```\n"
        )


def setup(app):
    app.connect("source-read", _local_toctree)
