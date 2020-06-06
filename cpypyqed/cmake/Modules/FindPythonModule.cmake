# Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
# From http://ompl.hg.sourceforge.net/hgweb/ompl/ompl/file/tip/CMakeModules/FindPython.cmake

##################################################
# Find python modules
##################################################
macro(find_python_module module name)
  set(module_name ${module})
  string(REGEX REPLACE \\. _ module ${module})
  string(TOUPPER ${module} module_upper)
  if(NOT PY_${module_upper})
    if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
      set(${module}_FIND_REQUIRED TRUE)
    endif()
    # A module's location is usually a directory, but for binary modules
    # it's a .so file.
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
                      "import re, ${module_name}; print re.compile('/__init__.py.*').sub('',${module_name}.__file__)"
                    RESULT_VARIABLE _${module}_status
                    OUTPUT_VARIABLE _${module}_location
                    ERROR_VARIABLE error
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _${module}_status)
      set(PY_${module_upper} ${_${module}_location} CACHE STRING
          "Location of Python module ${module}")
    endif(NOT _${module}_status)
  endif(NOT PY_${module_upper})
  find_package_handle_standard_args(${name} DEFAULT_MSG PY_${module_upper})
endmacro(find_python_module)
