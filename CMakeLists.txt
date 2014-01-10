cmake_minimum_required (VERSION 2.8.9)
project(cppqed)

enable_testing()

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

set(DOC_INSTALL_DIR ${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID})

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})

find_package(Boost QUIET COMPONENTS python)
if(Boost_PYTHON_FOUND)
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

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${DOC_INSTALL_DIR}
        PATTERN .git* EXCLUDE
        PATTERN .bzr EXCLUDE
        PATTERN *.kdev* EXCLUDE
        PATTERN build* EXCLUDE
)
add_subdirectory(Testing)

##################################################
# Documentation
##################################################

include(${PROJECT_SOURCE_DIR}/CPPQEDcore/cmake/Modules/CPPQEDUse.cmake)

cppqed_documentation("" core_doc elements_doc cpypyqed_doc)


##################################################
# Debian build
# if the directory debianbuild exists, configure the
# debian package build
##################################################

if(EXISTS ${PROJECT_SOURCE_DIR}/debianbuild)
  add_subdirectory(debianbuild)
endif()
