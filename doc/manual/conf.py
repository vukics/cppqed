# -*- coding: utf-8 -*-

execfile("../confCommon.py")

todo_include_todos=True

# The master toctree document.
master_doc = 'manual'

# General information about the project.
project = u'C++QEDv2 Reference Manual'


html_title = u'C++QED Reference Manual'

htmlhelp_basename = 'CQEDv2ReferenceManual'


latex_documents = [
  ('manual', 'C++QED_ReferenceManual.tex', u'C++QEDv2 Reference Manual',
   u'András Vukics', 'manual'),
  ('tutorial', 'C++QED_Tutorial.tex', ur'C++QEDv2 User Guide',
   ur'András Vukics', 'howto'),
  ('structureTutorial', 'C++QED_structureTutorial.tex', ur'C++QEDv2 Element Guide',
   ur'András Vukics', 'howto'),
]

