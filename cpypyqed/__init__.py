import sys

try:
  from core@DEBUG_SUFFIX@ import *
  from elements@DEBUG_SUFFIX@ import *
  from compilation.composite import *
except ImportError as e:
  ie = ImportError(e.message+ ". Maybe the wrong build type? This is cpypyqed configured for @CONF@.\n")
  raise ie
