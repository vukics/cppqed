cmake_minimum_required (VERSION 2.8.9)

project(cpypyqed)

include(FeatureSummary)

get_filename_component(CPPQED_CMAKE_MODULE_PATH "cmake/Modules" ABSOLUTE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CPPQED_CMAKE_MODULE_PATH})

find_package(CPPQED 2.99 REQUIRED)
include(${CPPQED_USE})
find_package(Boost QUIET REQUIRED COMPONENTS python)
find_package(PythonInterp 2 REQUIRED)
find_package(PythonLibs 2 REQUIRED)

include(PythonSetup)
CPPQED_SETUP()

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

SET(MANDATORY_LIBRARIES ${CPPQED_LIBRARIES} ${Boost_PYTHON_LIBRARY})

option(PYIO  "Python I/O" ON)
if(PYIO)
  set(ENABLE_PYTHON_IO On)
else(PYIO)
  set(ENABLE_PYTHON_IO Off)
endif(PYIO)

# Numpy
find_package(Numpy)
if( NOT NUMPY_FOUND )
  message(STATUS "Numpy not found, python I/O disabled.")
  set(ENABLE_PYTHON_IO Off)
endif( NOT NUMPY_FOUND )
# Check if python I/O module can be built, inform user
if ( ENABLE_PYTHON_IO )
  message(STATUS "Python input/output enabled.")
endif( ENABLE_PYTHON_IO )

set(PYTHON_MAX_RANK 8 CACHE STRING "Pre-instantiate template classes up to this rank in Python bindings.")
set_property(CACHE PYTHON_MAX_RANK PROPERTY STRINGS 2 4 6 8 10)

set(CPYPYQED_INSTALL_PATH_RELEASE ${PYTHON_MODULES_PATH}/cpypyqed)
set(CPYPYQED_INSTALL_PATH_DEBUG ${PYTHON_MODULES_PATH}/cpypyqed_d)

# Macro used for each module to select the correct extension according to the OS
MACRO(PYTHON_CPPMODULE arg1)
  SET_TARGET_PROPERTIES(${arg1} PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG   ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_d)
  IF(PYDEXTENSION)
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES SUFFIX ".pyd")
  ELSEIF(SOEXTENSION)
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES SUFFIX ".so")
  ENDIF(PYDEXTENSION)
  TARGET_LINK_LIBRARIES(${arg1} ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
  INSTALL(TARGETS ${arg1}
                RUNTIME DESTINATION ${CPYPYQED_INSTALL_PATH_RELEASE}
                LIBRARY DESTINATION ${CPYPYQED_INSTALL_PATH_RELEASE}
                CONFIGURATIONS Release
  )
  INSTALL(TARGETS ${arg1}
                RUNTIME DESTINATION ${CPYPYQED_INSTALL_PATH_DEBUG}
                LIBRARY DESTINATION ${CPYPYQED_INSTALL_PATH_DEBUG}
                CONFIGURATIONS Debug
  )
ENDMACRO(PYTHON_CPPMODULE)

# Macro used for sub-package creation
MACRO(PYTHON_SUBPACKAGE PACKAGE)
  SET(PACKAGENAME ${PACKAGE})
  AUX_SOURCE_DIRECTORY(. PY_SOURCE_FILES)

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

  configure_file(${PROJECT_SOURCE_DIR}/helpers/build_module.cc.in build${PACKAGENAME}.cc @ONLY)

  add_library(${PACKAGENAME} SHARED build${PACKAGENAME}.cc ${PY_SOURCE_FILES} ${ARGN})
  PYTHON_CPPMODULE(${PACKAGE})

ENDMACRO(PYTHON_SUBPACKAGE)

configure_file(${PROJECT_SOURCE_DIR}/helpers/config.h.in ${PROJECT_BINARY_DIR}/cpypyqed_config.h @ONLY)

include_directories(include)
include_directories(${PROJECT_BINARY_DIR})

# build python io module
if(ENABLE_PYTHON_IO)
  add_subdirectory(io)
endif()
add_subdirectory(core)
add_subdirectory(elements)

add_custom_target(cpypyqed DEPENDS io core elements)

set(INPUTS __init__.py compilation/cpypyqed_config.py compilation/ActTemplate.cc
           compilation/CompositeTemplate.cc compilation/composite.py compilation/ondemand.py
           compilation/__init__.py compilation/CMakeListsTemplate.txt)

foreach(file ${INPUTS})
  add_custom_command(TARGET cpypyqed PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${file}
      -DOUTPUT=${PROJECT_NAME}$<$<CONFIG:Debug>:_d>/${file}
      -DDEBUG_SUFFIX="$<$<CONFIG:Debug>:_d>"
      -DCONF=$<CONFIGURATION>
      -DCPPQED_CORE_GIT=${CPPQED_CORE_GIT}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
  )
  get_filename_component(DIR ${file} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}/${file}
          DESTINATION ${CPYPYQED_INSTALL_PATH_RELEASE}/${DIR}
          CONFIGURATIONS Release
  )
  install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_d/${file}
          DESTINATION ${CPYPYQED_INSTALL_PATH_DEBUG}/${DIR}
          CONFIGURATIONS Debug
  )
endforeach()


##################################################
# Documentation
##################################################

find_package(Doxygen)
find_package(Sphinx)
find_python_module(sphinxcontrib.doxylink)

if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND AND SPHINX_EXECUTABLE AND PY_SPHINXCONTRIB_DOXYLINK)
  set(CONF_DOC_DIR ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})
  file(MAKE_DIRECTORY ${CONF_DOC_DIR})
  file(RELATIVE_PATH CPPQED_RELATIVE_DOXYGEN_TAG ${CONF_DOC_DIR} ${CPPQED_DOXYGEN_TAG})
  get_filename_component(CPPQED_RELATIVE_DOXYGEN_DIR ${CPPQED_RELATIVE_DOXYGEN_TAG} PATH)
  set(CPPQED_RELATIVE_DOXYGEN_DIR ../${CPPQED_RELATIVE_DOXYGEN_DIR})
  configure_file(${PROJECT_SOURCE_DIR}/doc/conf.py.in ${CONF_DOC_DIR}/conf.py @ONLY)
  if(CPPQED_MONOLITHIC)
    configure_file(${PROJECT_SOURCE_DIR}/doc/layout.html.in ${CONF_DOC_DIR}/_templates/layout.html @ONLY)
  endif()
  foreach(AUXFILE index.rst userguide.rst)
    configure_file(${PROJECT_SOURCE_DIR}/doc/${AUXFILE} ${CONF_DOC_DIR} COPYONLY)
  endforeach()
  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
  # HTML output directory
  set(SPHINX_HTML_DIR "${CONF_DOC_DIR}/html")
  add_custom_target(${PROJECT_NAME}_doc
      ${SPHINX_EXECUTABLE}
          -q -b html
          -c "${CONF_DOC_DIR}"
          -d "${SPHINX_CACHE_DIR}"
          "${CONF_DOC_DIR}"
          "${SPHINX_HTML_DIR}"
      WORKING_DIRECTORY ${CONF_DOC_DIR}
      COMMENT "Building HTML documentation with Sphinx"
      DEPENDS doc/conf.py.in
  )
  add_custom_command(TARGET ${PROJECT_NAME}_doc PRE_BUILD
                     COMMAND ${CMAKE_COMMAND}
                        -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/doc/reference.rst.in
                        -DOUTPUT=${CONF_DOC_DIR}/reference.rst
                        -DMODULE_NAME=${PROJECT_NAME}$<$<CONFIG:Debug>:_d>
                        -P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
                     VERBATIM
  )
  add_dependencies(${PROJECT_NAME}_doc cpypyqed)
  if(CPPQED_MONOLITHIC)
    add_dependencies(${PROJECT_NAME}_doc core_doc elements_doc)
    install(DIRECTORY ${CONF_DOC_DIR}/html
          DESTINATION ${CPPQED_DOC_DIR}/${PROJECT_NAME}
          OPTIONAL
    )
  endif()
endif()


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
                                PURPOSE "Data input/output with Python.")

add_feature_info("Python input/output" ENABLE_PYTHON_IO "Python C++QED input/output module.")

feature_summary( WHAT ALL )
