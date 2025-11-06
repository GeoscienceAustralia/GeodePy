# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GeodePy'
copyright = '2025, Geoscience Australia'
author = 'Geoscience Australia'
release = '0.6.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_sidebars = {
    "**": [
        "sidebar/brand.html",      # Logo
        "sidebar/extra.html",      # Custom tagline + badge
        "sidebar/scroll-start.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/scroll-end.html",
    ]
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_theme_options = {
    "light_logo": "geodepy-logo-light.png",
    "dark_logo": "geodepy-logo-dark.png",
    "sidebar_hide_name": True,
}

html_css_files = [
    'custom.css',
]