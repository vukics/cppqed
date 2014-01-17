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

cmake_minimum_required (VERSION 2.8.9)
project(cppqed)

enable_testing()

include(GNUInstallDirs)

#! \name Project variables
#! These variables are used in the subprojects.
#! @{

#! \brief Other subprojects use this variable to determine if this is a monolithic build
set(CPPQED_MONOLITHIC 1)

#! \brief Install directory of the Doxygen documentation.
#!
#! Note that doxygen documentation can only be installed in monolithic builds. This
#! variable is used in cppqed_documentation().
set(CPPQED_DOC_DIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID}")

#! @}

#! \file
#! <!--#########################################################-->
#! ### Compilation of the components
#! <!--#########################################################-->
#!
#! Compile the components core, elements, scripts (optional) and cpypyqed (optional)
#! by adding their sub-directories.
#!
#! Also include the test suite.

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})

#! \name CMake options
#! These options can be switched with `-D<OPTION>=(ON|OFF)` when calling CMake.
#! @{

#! This CMake option determines if the Python modules should be compiled.
option(COMPILE_CPYPYQED "Compile Python wrapper and Python I/O module" On)
find_package(Boost QUIET COMPONENTS python)
if(Boost_PYTHON_FOUND AND COMPILE_CPYPYQED)
  add_subdirectory(cpypyqed)
endif()

#! \brief This CMake option determines if the scripts should be compiled as part
#!    of the monolithic build.
option(COMPILE_SCRIPTS "Compile the example scripts" On)
#! @}
if(COMPILE_SCRIPTS)
  add_subdirectory(CPPQEDscripts)
endif()

add_subdirectory(Testing)

#! \file
#! <!--#########################################################-->
#! ### Documentation
#! <!--#########################################################-->
#!
#! This section builds the Doxygen documentation (if Doxygen is found) and furthermore
#! installs the example projects "CustomElementsExample" and "CustomScriptsExample" to
#! the system.
#!
#! For the Doxygen documentation, the filter [CMakeDoxygenFilter](https://github.com/saschazelzer/CMakeDoxygenFilter),
#! is compiled, which allows to document %CMake in Doxygen. The overall C++QED documentation
#! is generated with cppqed_documentation(), and a target "doc" for the full documentation is created,
#! which depends on all the component documentation targets.

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  include(${core_SOURCE_DIR}/cmake/Modules/CPPQEDUse.cmake)
  add_executable(CMakeDoxygenFilter doc/CMakeDoxygenFilter.cpp)
  set_target_properties(CMakeDoxygenFilter PROPERTIES COMPILE_FLAGS -DUSE_NAMESPACE=CMake)
  set(CMakeDoxygenFilter_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/CMakeDoxygenFilter${CMAKE_EXECUTABLE_SUFFIX}")
  cppqed_documentation("cppqed_" "${CMAKE_BINARY_DIR}/doc/core/core.tag" core_doc CMakeDoxygenFilter)
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

feature_summary( WHAT ALL )
