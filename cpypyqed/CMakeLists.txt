cmake_minimum_required (VERSION 2.8.9)

project(cpypyqed)

# stuff should be built inside a "cpypyqed" subdirectory, so that we can
# easily import the build tree in python for testing
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

include(FeatureSummary)

get_filename_component(CPPQED_CMAKE_MODULE_PATH "cmake/Modules" ABSOLUTE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CPPQED_CMAKE_MODULE_PATH})

find_package(CPPQED 2.99 REQUIRED)
include(${CPPQED_USE})
find_package(Boost REQUIRED COMPONENTS python)
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
else()
  include_directories(SYSTEM ${NUMPY_INCLUDE_DIRS})
endif( NOT NUMPY_FOUND )
# Check if python I/O module can be built, inform user
if ( ENABLE_PYTHON_IO )
  message(STATUS "Python input/output enabled.")
endif( ENABLE_PYTHON_IO )

set(CPYPYQED_INSTALL_PATH ${PYTHON_MODULES_PATH}/cpypyqed)

set(PYTHON_MAX_RANK 11 CACHE STRING "Pre-instantiate template classes up to this rank in Python bindings.")
set_property(CACHE PYTHON_MAX_RANK PROPERTY STRINGS 1 2 3 4 5 6 7 8 9 10 11 12)

# Macro used for each module to select the correct extension according to the OS
MACRO(PYTHON_CPPMODULE arg1)
  SET_TARGET_PROPERTIES(${arg1} PROPERTIES PREFIX "" DEBUG_POSTFIX "")
  IF(PYDEXTENSION)
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES SUFFIX ".pyd")
  ELSEIF(SOEXTENSION)
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES SUFFIX ".so")
  ENDIF(PYDEXTENSION)
  TARGET_LINK_LIBRARIES(${arg1} ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
  INSTALL(TARGETS ${arg1}
                RUNTIME DESTINATION ${CPYPYQED_INSTALL_PATH}
                LIBRARY DESTINATION ${CPYPYQED_INSTALL_PATH}
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

add_custom_target(cpypyqed)

# build python io module
if(ENABLE_PYTHON_IO)
  add_subdirectory(io)
endif()
add_subdirectory(core)

foreach(file __init__.py)
  configure_file(${file} ${LIBRARY_OUTPUT_PATH}/${file} COPYONLY)
  install(FILES ${file} DESTINATION ${CPYPYQED_INSTALL_PATH})
endforeach(file)

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
