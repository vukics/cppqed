# -*- coding: utf-8 -*-

execfile("confCommon.py")

todo_include_todos=False

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'C++QED'


html_title = u'C++QED Documentation'

htmlhelp_basename = 'CQEDdoc'


html_theme_options = {
    "footerbgcolor" : "#808080",
    "footertextcolor" : "black",
    "sidebarbgcolor" : "#fff0f0",
    "sidebartextcolor" : "#800000",
    "sidebarlinkcolor" : "#804040",
    "relbarbgcolor" : "#800000",
    "relbartextcolor" : "#804040",
    "relbarlinkcolor" : "#fff0f0",
    "linkcolor" : "#804040",
    "visitedlinkcolor" : "#804040",
    "headbgcolor" : "#f0f0f0",
    "headtextcolor" : "#404040",
    "headlinkcolor" : "#808080",
    "bodyfont" : "Palatino Linotype",
    "headfont" : "Palatino Linotype"
}

latex_documents = [
  ('index', 'C++QED.tex', u'C++QED Documentation',
   u'András Vukics', 'manual'),
]

latex_use_parts = True
