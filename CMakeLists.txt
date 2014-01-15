cmake_minimum_required (VERSION 2.8.9)
project(cppqed)

enable_testing()

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

set(CPPQED_DOC_DIR ${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID})

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})

option(COMPILE_CPYPYQED "Compile Python wrapper and Python I/O module" On)
find_package(Boost QUIET COMPONENTS python)
if(Boost_PYTHON_FOUND AND COMPILE_CPYPYQED)
  add_subdirectory(cpypyqed)
endif()

# Disabled CustomExample compilation by default for now, add it as a test later
# add_subdirectory(CustomElementsExample)
# set(CPPQEDelements_custom_example_DIR ${elements_custom_example_BINARY_DIR})
# add_subdirectory(CustomScriptsExample)

option(COMPILE_SCRIPTS "Compile the example scripts" On)
if(COMPILE_SCRIPTS)
  add_subdirectory(CPPQEDscripts)
endif()

add_subdirectory(Testing)

##################################################
# Documentation
##################################################

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  include(${core_SOURCE_DIR}/cmake/Modules/CPPQEDUse.cmake)
  add_executable(CMakeDoxygenFilter doc/CMakeDoxygenFilter.cpp)
  set_target_properties(CMakeDoxygenFilter PROPERTIES COMPILE_FLAGS -DUSE_NAMESPACE=CMake)
  set(CMakeDoxygenFilter_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/CMakeDoxygenFilter${CMAKE_EXECUTABLE_SUFFIX}")
  cppqed_documentation("cppqed_" "${CMAKE_BINARY_DIR}/doc/core/core.tag" core_doc)
  add_dependencies(cppqed_doc CMakeDoxygenFilter)
  add_custom_target(doc)
  add_dependencies(doc cppqed_doc core_doc elements_doc cpypyqed_doc)
endif()

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${CPPQED_DOC_DIR}
        PATTERN .git* EXCLUDE
        PATTERN .bzr EXCLUDE
        PATTERN *.kdev* EXCLUDE
        PATTERN build* EXCLUDE
)

##################################################
# Debian build
# if the directory debianbuild exists, configure the
# debian package build
##################################################

if(EXISTS ${PROJECT_SOURCE_DIR}/debianbuild)
  add_subdirectory(debianbuild)
endif()

feature_summary( WHAT ALL )
