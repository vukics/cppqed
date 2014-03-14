import sys
import cpypyqed.config

from ..core_d import *
from ..elements_d import *
cpypyqed.config.build_type="debug"
cpypyqed.config.module_suffix="_d"
from ..compilation.composite import *
sys.stderr.write("Warning: Using C++QED libraries with Debug configuration.\n")
sys.stderr.write("Warning: Simulations will be slow.\n")

