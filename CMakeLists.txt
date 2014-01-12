cmake_minimum_required (VERSION 2.8.9)
project(cppqed)

enable_testing()

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

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

include(${core_SOURCE_DIR}/cmake/Modules/CPPQEDUse.cmake)
cppqed_documentation("cppqed_" "${CMAKE_BINARY_DIR}/doc/core/core.tag" core_doc)

add_custom_target(doc)
add_dependencies(doc cppqed_doc core_doc elements_doc cpypyqed_doc)

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${DOC_INSTALL_DIR}
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
