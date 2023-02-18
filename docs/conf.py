import subprocess, os

read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"
build_directory = "build"
if read_the_docs_build:
    subprocess.run(
        ["doxygen"], shell=True, env={"DOXYGEN_OUTPUT_DIRECTORY": build_directory}
    )

# -- Project information -----------------------------------------------------

project = "Catima"
copyright = "2017, A. Prochazka"
author = "A. Prochazka"


# -- General configuration ---------------------------------------------------

extensions = ["recommonmark", "breathe"]

breathe_projects = {"catima": os.path.join(build_directory, "doxygenxml")}
breathe_default_project = "catima"
breathe_domain_by_extenstion = {"h": "cpp"}
primary_domain = "cpp"
highlight_language = "cpp"
# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
