# Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \addtogroup Main
#! @{

#! \file
#! \brief Top level %CMake file controlling monolithic builds.
#!
#! In monolithic builds, the paths to the subprojects (`CPPQEDcore`, `CPPQEDelements`, `CPPQEDscripts`, `cpypyqed`)
#! are known. All these components are built with calls to `add_subdirectory` from within this file. The subprojects
#! are built just as if they were standalone projects (with very few exceptions where the variable `CPPQED_MONOLITHIC`
#! is queried and things are handled differently if it is defined). `CPPQEDcore` and `CPPQEDelements`
#! export their relevant targets and other subprojects import them, just as for standalone projects. However, this
#! file sets the variables `CPPQED_DIR` and `CPPQEDelements_DIR`
#! (c.f. \ref cmake_find_components "how CMake finds components") to the build directories of this monolithic build,
#! so that the right subprojects are found even if C++QED is installed or other build directories are registered.
#!
#! This CMakeLists file has the following structure:

#! @}

cmake_minimum_required (VERSION 3.12)
project(cppqed)

enable_testing()

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_SOURCE_DIR}/CPPQEDcore/cmake/Modules" "${CMAKE_SOURCE_DIR}/cpypyqed/cmake/Modules")

#! \name Project variables
#! These variables are used in the subprojects.
#! @{

#! \brief Other subprojects use this variable to determine if this is a monolithic build
set(CPPQED_MONOLITHIC 1)

set(cppqed_DOC_DIR ${CMAKE_BINARY_DIR}/doc/cppqed)
set(core_DOC_DIR ${CMAKE_BINARY_DIR}/doc/core)
set(core_DOXYGEN_TAG ${core_DOC_DIR}/core.tag)
set(elements_DOC_DIR ${CMAKE_BINARY_DIR}/doc/elements)
set(elements_DOXYGEN_TAG ${elements_DOC_DIR}/elements.tag)
set(cpypyqed_DOC_DIR ${CMAKE_BINARY_DIR}/doc/cpypyqed)

#! @}

find_package(Boost QUIET COMPONENTS python)
find_package(PythonInterp 3 QUIET)
if(${PYTHON_VERSION_STRING} VERSION_LESS 3.5)
  message(WARNING "Python >=3.5 needed, disabling cpypyqed and testsuite.")
  set(PYTHONINTERP_FOUND 0)
else()
  find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT QUIET)
endif()
find_package(Numpy QUIET)


#! \file
#! <!--#########################################################-->
#! ### Compilation of the components
#! <!--#########################################################-->
#!
#! Compile the components core, elements, scripts (optional) and cpypyqed (optional)
#! by adding their sub-directories.

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})

#! \brief Install directory of the Doxygen documentation.
#!
#! Note that doxygen documentation can only be built and installed in monolithic builds.
set(CPPQED_DOC_DIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID}")

#! \name Project options
#! @{

#! This %CMake option determines if the Python modules should be compiled.
option(COMPILE_CPYPYQED "Compile Python wrapper and Python I/O module" On)

# For some reason kdevelop does not recognize PYTHONLIBS_FOUND, use PYTHON_INCLUDE_DIRS instead
if(Boost_PYTHON_FOUND AND PYTHONINTERP_FOUND AND PYTHON_INCLUDE_DIRS AND NUMPY_FOUND AND COMPILE_CPYPYQED)
  set(ENABLE_CPYPYQED 1)
endif()
add_subdirectory(cpypyqed)


#! \brief This %CMake option determines if the scripts should be compiled as part
#!    of the monolithic build.
option(COMPILE_SCRIPTS "Compile the example scripts" On)

#! \brief This %CMake configuration variable lets users exclude scripts from the ALL target. These
#!    scripts can still be compiled by directly calling their target.
set(CPPQED_EXCLUDE_SCRIPTS "" CACHE STRING "Exclude from the ALL target")
#! @}
if(COMPILE_SCRIPTS)
  add_subdirectory(CPPQEDscripts)
endif()

add_subdirectory(Testing)

# The EXCLUDE_FROM_ALL target makes Xcode scip the subdirectory entirely
if(XCODE)
  add_subdirectory(CustomElementsExample)
  set(CPPQEDelements_custom_example_DIR ${elements_custom_example_BINARY_DIR})
  add_subdirectory(CustomScriptsExample)
else()
  add_subdirectory(CustomElementsExample EXCLUDE_FROM_ALL)
  set(CPPQEDelements_custom_example_DIR ${elements_custom_example_BINARY_DIR})
  add_subdirectory(CustomScriptsExample EXCLUDE_FROM_ALL)
endif()

#! \file
#! <!--#########################################################-->
#! ### Source package
#! <!--#########################################################-->
#!
#! Create source tar ball with CPack for publication on sourceforge.

set(CPACK_PACKAGE_VERSION_MAJOR ${CPPQED_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${CPPQED_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${CPPQED_VERSION_PATCH})

set(CPACK_GENERATOR TGZ)

set(CPACK_SOURCE_PACKAGE_FILE_NAME "cppqed-${CPPQED_VERSION}" CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_IGNORE_FILES "~$" "/build/" "\\\\.kdev4.*" "\\\\.git.*" ".*orig$")

include(CPack)

#! \file
#! <!--#########################################################-->
#! ### Documentation
#! <!--#########################################################-->
#!
#! This section builds the Doxygen documentation (if Doxygen is found) and furthermore
#! installs the example projects "CustomElementsExample" and "CustomScriptsExample" to
#! the system. The install location is \ref CMake::CPPQED_DOC_DIR "CPPQED_DOC_DIR".
#!
#! For the Doxygen documentation, the filter [CMakeDoxygenFilter](https://github.com/saschazelzer/CMakeDoxygenFilter),
#! is compiled, which allows to document %CMake in Doxygen. The overall C++QED documentation
#! is generated with cppqed_documentation(), and a target "doc" for the full documentation is created,
#! which depends on all the component documentation targets.

include(${core_SOURCE_DIR}/cmake/Modules/CPPQEDUse.cmake)
add_executable(CMakeDoxygenFilter doc/CMakeDoxygenFilter.cpp)
set_target_properties(CMakeDoxygenFilter PROPERTIES COMPILE_FLAGS -DUSE_NAMESPACE=CMake)
set(CMakeDoxygenFilter_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/CMakeDoxygenFilter${CMAKE_EXECUTABLE_SUFFIX}")
cppqed_documentation("cppqed_" "${CMAKE_BINARY_DIR}/doc/core/core.tag;${CMAKE_BINARY_DIR}/doc/elements/elements.tag" core_doc elements_doc CMakeDoxygenFilter)
if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  add_custom_target(doc)
  add_dependencies(doc cppqed_doc core_doc elements_doc)
  if(DOXYLINK_FOUND)
    add_dependencies(cpypyqed_doc)
  endif()
  configure_file(doc/index.html ${CMAKE_BINARY_DIR}/doc COPYONLY)
endif()

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${CPPQED_DOC_DIR}
        PATTERN .git* EXCLUDE
        PATTERN .bzr EXCLUDE
        PATTERN *.kdev* EXCLUDE
        PATTERN build* EXCLUDE
)

#! \file
#! <!--#########################################################-->
#! ### Debian build
#! <!--#########################################################-->
#!
#! This adds the sub-directory "debianbuild" if it exists. If not this section is silently ignored.
#! The directory "debianbuild" is not part of the C++QED project. One can put a special
#! [repository](https://ge-c705.uibk.ac.at/cppqed/debian) here, which will then generate debian
#! directories for building source packages. This guarantees that debian packages have the right
#! version numbers and soversion in their package names, and also helps to maintain packages for different
#! distributions without code duplication.

if(EXISTS ${PROJECT_SOURCE_DIR}/debianbuild)
  add_subdirectory(debianbuild)
endif()

add_feature_info("cpypyqed" ENABLE_CPYPYQED "Python part of C++QED (wrapper and input/output).")
feature_summary( WHAT ALL )
