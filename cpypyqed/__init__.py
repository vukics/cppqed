import sys
import cpypyqed.config

if '--debug' in sys.argv:
  from .debug import *
  sys.argv.remove('--debug')
else:
  try:
    from core import *
    from elements import *
    cpypyqed.config.build_type="release"
    cpypyqed.config.module_suffix=""
    from compilation.composite import *

  except ImportError as e:
    from .debug import *
