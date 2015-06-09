# Copyright Raimar Sandner 2012–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \ingroup Main
#! \file
#! \brief Top level %CMake file for the C++QED cpypyqed component.
#!
#! The file has the following structure:


cmake_minimum_required (VERSION 2.8.9)

project(cpypyqed)

include(FeatureSummary)

if(CPPQED_MONOLITHIC AND NOT ENABLE_CPYPYQED)
  return()
endif()

#! \file
#! <!--#########################################################-->
#! ### Dependencies
#! <!--#########################################################-->
#!
#! Find dependencies for cpypyqed, i.e. C++QED, boost-python, the python interpreter,
#! python libraries and numpy.

get_filename_component(CPPQED_CMAKE_MODULE_PATH "cmake/Modules" ABSOLUTE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CPPQED_CMAKE_MODULE_PATH})
message(STATUS "Debug: ${PYTHON_EXECUTABLE}")

find_package(CPPQED 2.100 REQUIRED)
include(${CPPQED_USE})
if(NOT PYTHONINTERP_FOUND)
  find_package(PythonInterp 2.6 REQUIRED)
endif()
find_package(PythonLibs 2.6 REQUIRED)
find_package(Boost QUIET REQUIRED COMPONENTS python)
find_package(Numpy REQUIRED)

include(PythonSetup)
CPPQED_SETUP()
include_directories(${CPPQED_INCLUDE_DIRS})

include_directories(${CPPQED_INCLUDE_DIRS})
include_directories(SYSTEM ${NUMPY_INCLUDE_DIRS})
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

SET(MANDATORY_LIBRARIES ${CPPQED_LIBRARIES} ${Boost_PYTHON_LIBRARY})

#! \file
#! <!--#########################################################-->
#! ### Configuration
#! <!--#########################################################-->
#!
#! Set \ref CMake::PYTHON_MAX_RANK "PYTHON_MAX_RANK" and the install paths where the module should end up.
#! Note that in debug build mode a `_d` suffix will be appended to the
#! module name. This way debug and release modules can be installed
#! simultaneously.
#!
#! The \ref CMake::PYTHON_MAX_RANK "PYTHON_MAX_RANK" variable ends up in a generated `cpypyqed_config.h`
#! header file and can be included where needed.

#! \name Project options
#! @{

#! \brief The maximum rank up to which some class templates with RANK parameter should be
#!   pre-initialized (must be even).
set(PYTHON_MAX_RANK 8 CACHE STRING "Pre-instantiate template classes up to this rank in Python bindings.")
#! @}
set_property(CACHE PYTHON_MAX_RANK PROPERTY STRINGS 2 4 6 8 10)

set(CPYPYQED_INSTALL_PATH ${PYTHON_MODULES_PATH}/cpypyqed)

#! \file
#! <!--#########################################################-->
#! ### Helper macros
#! <!--#########################################################-->
#!
#! The helper macros python_cppmodule() and python_submodule().

#! \ingroup Helpers
#! \brief This macro selects the correct extension for Python modules according to the OS.
#! \param modulename The name of the Python module target.
#!
#! It also links the module to whatever is currently set in the variable `MANDATORY_LIBRARIES`
#! and the Python libraries, and installs the module.
macro(python_cppmodule modulename)
  SET_TARGET_PROPERTIES(${modulename} PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PROJECT_NAME}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
  IF(PYDEXTENSION)
    SET_TARGET_PROPERTIES(${modulename} PROPERTIES SUFFIX ".pyd")
  ELSEIF(SOEXTENSION)
    SET_TARGET_PROPERTIES(${modulename} PROPERTIES SUFFIX ".so")
  ENDIF(PYDEXTENSION)
  TARGET_LINK_LIBRARIES(${modulename} ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
  INSTALL(TARGETS ${modulename}
                RUNTIME DESTINATION ${CPYPYQED_INSTALL_PATH}
                LIBRARY DESTINATION ${CPYPYQED_INSTALL_PATH}
  )
endmacro()

#! \ingroup Helpers
#! \brief Macro used for Python submodule creation.
#! \param modulename The name of the Python module (without extension)
#!
#! This macro is meant to be called from a `CMakeLists.txt` file residing in sub-directories of
#! the top-level cpypyqed directory.
#!
#! This macro will take all source files in the current directory and search for functions
#! matching `export_.*`. The names of these functions are written into a boost preprocessor
#! sequence and saved in the variable `EXPORT_HELPER_SEQUENCE` (sorted alphabetically).
#!
#! From the template file `helpers/build_module.cc.in`, a source file is generated, in which this
#! preprocessor sequence is used to declare and call all the export functions. A target is created
#! building the generated `build${modulename}.cc` together with all source files in the current
#! directory, and linking to the Python and boost libraries. Additional libraries needed
#! (e.g. elements libraries) have to be added to the variable `MANDATORY_LIBRARIES` before calling
#! this macro.
#!
#! Finally python_cppmodule() is called for correct linkage and installation.
macro(python_submodule modulename)
  set(modulename ${modulename})
  aux_source_directory(. PY_SOURCE_FILES)

  foreach(s ${PY_SOURCE_FILES})
    file( STRINGS ${s} EXPORT_HELPER_LINE REGEX "[ ]*void[ ]+export_" )
    string(REGEX REPLACE ".*(export_.*)\\(.*" "\\1" EXPORT_HELPER "${EXPORT_HELPER_LINE}"  )
    set(EXPORT_HELPERS ${EXPORT_HELPERS} ${EXPORT_HELPER})
  endforeach(s)
  if(EXPORT_HELPERS)
    list(SORT EXPORT_HELPERS)
  endif(EXPORT_HELPERS)
  foreach(s ${EXPORT_HELPERS})
    set(EXPORT_HELPER_SEQUENCE ${EXPORT_HELPER_SEQUENCE}\(${s}\))
  endforeach(s)

  configure_file(${PROJECT_SOURCE_DIR}/helpers/exportfunctions.h.in exportfunctions.h @ONLY)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  configure_file(${PROJECT_SOURCE_DIR}/helpers/build_module.cc.in build${modulename}.cc @ONLY)

  add_library(${modulename} SHARED build${modulename}.cc ${PY_SOURCE_FILES} ${ARGN})
  add_dependencies(${modulename} cpypyqed_setup)
  python_cppmodule(${modulename})

endmacro()

#! \file
#! <!--#########################################################-->
#! ### Compilation
#! <!--#########################################################-->
#!
#! Add the sub-directories io, core and elements, which basically only call python_submodule().

include_directories(include)
include_directories(${PROJECT_BINARY_DIR})

add_custom_target(cpypyqed_setup)

add_subdirectory(io)
add_subdirectory(core)
add_subdirectory(elements)

add_custom_target(cpypyqed ALL)
add_dependencies(cpypyqed io core elements)

add_custom_command(TARGET cpypyqed_setup PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DINPUT=${PROJECT_SOURCE_DIR}/helpers/config.h.in
      -DOUTPUT=cpypyqed_config.h
      -DDEBUG_SUFFIX="$<$<CONFIG:Debug>:_d>"
      -DPACKAGE_NAME="cpypyqed$<$<CONFIG:Debug>:.debug>"
      -DPYTHON_MAX_RANK=${PYTHON_MAX_RANK}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
  )

#! \file
#! <!--#########################################################-->
#! ### Prepare Python package
#! <!--#########################################################-->
#!
#! Copy all the files which are needed for the final Python package
#! to the build directory. The package name is `cpypyqed`. Because it is possible that the
#! build configuration is only available at build time (e.g. Mac OS X, XCode),
#! not at configure time, we cannot use `configure_file` directly.
#!
#! Instead we use a little trick: The call to `configure_file` is wrapped in a
#! trivial %CMake script `generate.cmake`, which is called at build time.
#! For this, `add_custom_command` is used to set up a pre-dependency for the
#! `cpypyqed` target. All variables which should be available for substitution
#! have to be passed in by `add_custom_command` with `-D` defines. In the custom
#! command we can make use of the %CMake generator expression, e.g. `$<$<CONFIG:Debug>:_d>`,
#! which expands to `_d` in Debug mode and does nothing in Release mode.
#!
#! This section also takes care to install the additional package files to the system.

set(INPUTS "__init__.py" "config.py" "compilation/ActTemplate.cc"
           "compilation/CompositeTemplate.cc" "compilation/composite.py" "compilation/ondemand.py"
           "compilation/__init__.py" "compilation/CMakeListsTemplate.txt" "debug/__init__.py" "tools/__init__.py" 
           "tools/quantumstate.py" "tools/expvalues.py" "tools/utils.py" "tools/initialconditions.py" "tools/io.py"
           "tools/test_statevector.py" "tools/passthrough.py" "tools/operators.py")

foreach(file ${INPUTS})
  add_custom_command(TARGET cpypyqed PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${file}
      -DOUTPUT=${PROJECT_NAME}/${file}
      -DDEBUG_SUFFIX="$<$<CONFIG:Debug>:_d>"
      -DCONF=$<CONFIGURATION>
      -DCPPQED_CORE_GIT=${CPPQED_CORE_GIT}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
  )
  get_filename_component(DIR ${file} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}/${file}
          DESTINATION ${CPYPYQED_INSTALL_PATH}/${DIR}
  )
endforeach()

#! \file
#! <!--#########################################################-->
#! ### Documentation
#! <!--#########################################################-->
#!
#! Find Doxygen, Sphinx and Doxylink. If all these components are found,
#! build the documentation. Sphinx has the advantage that it actually loads
#! the modules with a python interpreter and looks at the docstrings.
#! Otherwise it would be very cumbersome to document a boost python module,
#! because all the signatures would have to be duplicated in the documentation
#! code.
#!
#! The relative path to the Doxygen tagfile `core.tag` is established, in order
#! to link to the C++QED API with Doxylink.
#!
#! In monolithic mode, the file `doc`layout.html.in` is copied over to the
#! documentation build directory. This files contains links to the other C++QED
#! components and makes them appear on the Sphinx pages.
#!
#! The sphinx source file `reference.rst` depends on the Python package name,
#! therefore we again have to generate this file at build time, when the
#! build configuration is known in any case. The solution is again, as in the
#! previous section, to let a custom command generate the file with a simple
#! cmake script.
#!
#! Finally, the call to the sphinx binary is added as a custom command, the
#! target is named `cpypyqed_doc`. In monolithic mode, the resulting html
#! directory is installed to the system.

find_package(Doxygen 1.8)
find_package(Sphinx)
find_package(Doxylink)

if(CPPQED_MONOLITHIC AND DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND AND SPHINX_FOUND AND DOXYLINK_FOUND)
  set(CPPQED_DOC_DIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID}")
  file(MAKE_DIRECTORY ${${PROJECT_NAME}_DOC_DIR})
  file(RELATIVE_PATH core_RELATIVE_DOXYGEN_TAG ${${PROJECT_NAME}_DOC_DIR} ${core_DOXYGEN_TAG})
  get_filename_component(core_RELATIVE_DOXYGEN_DIR ${core_RELATIVE_DOXYGEN_TAG} PATH)
  configure_file(${PROJECT_SOURCE_DIR}/doc/conf.py.in ${${PROJECT_NAME}_DOC_DIR}/conf.py @ONLY)
  if(CPPQED_MONOLITHIC)
    configure_file(${PROJECT_SOURCE_DIR}/doc/layout.html.in ${${PROJECT_NAME}_DOC_DIR}/_templates/layout.html @ONLY)
  endif()
  foreach(AUXFILE index.rst userguide.rst tools.rst)
    configure_file(${PROJECT_SOURCE_DIR}/doc/${AUXFILE} ${${PROJECT_NAME}_DOC_DIR} COPYONLY)
  endforeach()
  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
  # HTML output directory
  set(SPHINX_HTML_DIR "${${PROJECT_NAME}_DOC_DIR}/html")
  add_custom_command(OUTPUT ${${PROJECT_NAME}_DOC_DIR}/reference.rst
                     COMMAND ${CMAKE_COMMAND}
                        -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/doc/reference.rst.in
                        -DOUTPUT=${${PROJECT_NAME}_DOC_DIR}/reference.rst
                        -DDEBUG_SUFFIX="$<$<CONFIG:Debug>:_d>"
                        -P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
                     DEPENDS generate.cmake doc/reference.rst.in
  )
  add_custom_target(${PROJECT_NAME}_doc
      ${SPHINX_EXECUTABLE}
          -q -b html
          -c "${${PROJECT_NAME}_DOC_DIR}"
          -d "${SPHINX_CACHE_DIR}"
          "${${PROJECT_NAME}_DOC_DIR}"
          "${SPHINX_HTML_DIR}"
      WORKING_DIRECTORY ${${PROJECT_NAME}_DOC_DIR}
      COMMENT "Building HTML documentation with Sphinx"
      DEPENDS ${${PROJECT_NAME}_DOC_DIR}/reference.rst
  )
  add_dependencies(${PROJECT_NAME}_doc cpypyqed)
  if(CPPQED_MONOLITHIC)
    add_dependencies(${PROJECT_NAME}_doc core_doc elements_doc)
    install(DIRECTORY ${${PROJECT_NAME}_DOC_DIR}/html
          DESTINATION ${CPPQED_DOC_DIR}/${PROJECT_NAME}
          OPTIONAL
    )
  endif()
endif()

#! \file
#! <!--#########################################################-->
#! ### Feature summary
#! <!--#########################################################-->
#!
#! Generate a summary of the cpypyqed dependencies.

set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                DESCRIPTION "Collection of portable C++ source libraries."
                                TYPE REQUIRED
                                PURPOSE "Python wrapper library for C++ libraries.")
set_package_properties(PythonLibs PROPERTIES URL "http://www.python.org/"
                                DESCRIPTION "General-purpose, high-level programming language."
                                TYPE OPTIONAL
                                PURPOSE "Access C++QED functionality through Python (work in progress, incomplete).")
set_package_properties(Numpy PROPERTIES URL "http://www.numpy.org/"
                                DESCRIPTION "NumPy is the fundamental package for scientific computing with Python"
                                TYPE OPTIONAL
                                PURPOSE "Data input/output with Python, testsuite.")
set_package_properties(Scipy PROPERTIES URL "http://www.scipy.org/"
                                DESCRIPTION "SciPy is a Python-based ecosystem of open-source software for mathematics, science, and engineering."
                                TYPE OPTIONAL
                                PURPOSE "The testsuite depends on this.")
set_package_properties(Sphinx PROPERTIES URL "http://sphinx-doc.org/"
                                DESCRIPTION "Sphinx is a documentation generator which converts reStructuredText files into HTML websites"
                                TYPE OPTIONAL
                                PURPOSE "Sphinx is needed to build the cpypyqed documentation.")
set_package_properties(Doxylink PROPERTIES URL "https://pypi.python.org/pypi/sphinxcontrib-doxylink"
                                DESCRIPTION "Doxylink is a Sphinx extension to link to external Doxygen API documentation."
                                PURPOSE "Needed to link to the rest of the C++QED in the cpypyqed documentation.")

if(NOT CPPQED_MONOLITHIC)
  feature_summary( WHAT ALL )
endif()
