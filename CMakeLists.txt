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

# Canopy workaround: the python interpreter is found, but not the libs and include dirs
if(PYTHONINTERP_FOUND AND NOT PYTHONLIBS_FOUND)
  execute_process(COMMAND "${PYTHON_EXECUTABLE}-config" "--includes"
    OUTPUT_VARIABLE _PYTHON_CONFIG_INCLUDES
    RESULT_VARIABLE _PYTHON_CONFIG_SUCCESS)
  if(_PYTHON_CONFIG_SUCCESS MATCHES 0)
    string(REGEX REPLACE "-I(.*) .*" "\\1" PYTHON_INCLUDE_DIR ${_PYTHON_CONFIG_INCLUDES})
    get_filename_component(_PYTHON_BASE ${PYTHON_INCLUDE_DIR} DIRECTORY)
    get_filename_component(_PYTHON_BASE ${_PYTHON_BASE} DIRECTORY)
    set(PYTHON_LIBRARY ${_PYTHON_BASE}/lib/libpython${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${CMAKE_SHARED_LIBRARY_SUFFIX})
    message(STATUS "Trying to find python library at ${PYTHON_LIBRARY}")
    message(STATUS "Trying to find python includes at ${PYTHON_INCLUDE_DIR}")
    find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT QUIET)
  endif()
endif()


#! \file
#! <!--#########################################################-->
#! ### Compilation of dependencies
#! <!--#########################################################-->
#!
#! This is a convenience feature to download and compile some selected
#! dependencies automatically. It uses %CMake's
#! [ExternalProject](http://www.cmake.org/cmake/help/v2.8.12/cmake.html#module:ExternalProject).
#!
#! We use a trick to install the bundled libraries. If the bundled version of a library is used,
#! we link it privately (\ref CMake::PRIVATE_LIBS instead of \ref CMake::PUBLIC_LIBS). Then
#! we add the library manually to `CPPQED_LIBRARIES` in CPPQEDConfig.cmake.in, with the correct
#! location for the build tree and installed tree, respectively.

#! If this option is set, %CMake downloads and compiles the blitz++ dependency automatically.
option(BUNDLED_BLITZ "Download and compile blitz++ automatically" Off)

if(BUNDLED_BLITZ)
  include(ExternalProject)

  find_package(Boost QUIET COMPONENTS serialization)
  if(NOT Boost_LIBRARY_DIR)
    set(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR_RELEASE})
  endif()
  if(Boost_SERIALIZATION_FOUND)
    set(BLITZ_CONFIG_OPTIONS "${BLITZ_CONFIG_OPTIONS} --enable-serialization ")
    set(blitz_SERIALIZATION_FOUND On)
  endif()

  string(REGEX REPLACE "/include" "" local_boost_root ${Boost_INCLUDE_DIR})

  if( CMAKE_SIZEOF_VOID_P EQUAL 8 ) # 64bit system
    set(BLITZ_CONFIG_OPTIONS "${BLITZ_CONFIG_OPTIONS} --enable-64bit ")
  else( CMAKE_SIZEOF_VOID_P EQUAL 8 ) # 32bit system
    set(BLITZ_CONFIG_OPTIONS "${BLITZ_CONFIG_OPTIONS} --enable-simd-width=8 ")
  endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

  find_program(autoreconf autoreconf)
  if(NOT autoreconf)
    message(FATAL_ERROR "Program autoreconf not found, needed to compile blitz.")
  endif()
  get_filename_component(AUTOTOOLS_PATH ${autoreconf} DIRECTORY)

  if(XCODE)
    set(AUTORECONF_NO_OUTPUT "&>/dev/null")
  endif()
  configure_file(misc/blitz-make.sh.in ${PROJECT_BINARY_DIR}/blitz-make.sh @ONLY)

  if(EXISTS ${CMAKE_SOURCE_DIR}/external/blitz)
    set(BLITZ_SOURCE DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/blitz/. <SOURCE_DIR>/)
  else()
    set(BLITZ_SOURCE HG_REPOSITORY http://hg.code.sf.net/p/cppqed/blitz)
  endif()

  ExternalProject_Add(
    cppqed-blitz
    ${BLITZ_SOURCE}
    UPDATE_COMMAND ""
    PATCH_COMMAND patch --dry-run -N -p1 -d <SOURCE_DIR> -i ${CMAKE_SOURCE_DIR}/misc/blitz_rename_library.patch | grep "previously applied" || patch -p1 -d <SOURCE_DIR> -i ${CMAKE_SOURCE_DIR}/misc/blitz_rename_library.patch
          COMMAND patch --dry-run -N -p1 -d <SOURCE_DIR> -i ${CMAKE_SOURCE_DIR}/misc/blitz_python26.patch | grep "previously applied" || patch -p1 -d <SOURCE_DIR> -i ${CMAKE_SOURCE_DIR}/misc/blitz_python26.patch
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/blitz-make.sh <SOURCE_DIR> <INSTALL_DIR>
    BUILD_COMMAND make lib
  )
  ExternalProject_Get_Property(cppqed-blitz install_dir)
  set(blitz_INCLUDE_DIRS ${install_dir}/include)
  set(blitz_LIBRARIES ${install_dir}/lib/libblitzcppqed${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(BUNDLED_BLITZ_INSTALLED_LIBRARY ${CMAKE_INSTALL_FULL_LIBDIR}/libblitzcppqed${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(blitz_FOUND 1)
endif()

#! If this option is set, %CMake downloads and compiles the optional FLENS dependency automatically.
option(BUNDLED_FLENS "Download and compile FLENS automatically" Off)

if(BUNDLED_FLENS)
  include(ExternalProject)

  if(EXISTS ${CMAKE_SOURCE_DIR}/external/flens)
    set(FLENS_SOURCE DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/flens/. <SOURCE_DIR>/)
  else()
    set(FLENS_SOURCE GIT_REPOSITORY git://github.com/michael-lehn/FLENS.git GIT_TAG "public")
  endif()

  ExternalProject_Add(
    flens
    ${FLENS_SOURCE}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ../flens/
    BUILD_COMMAND
    INSTALL_COMMAND make install
  )
  ExternalProject_Get_Property(flens install_dir)
  set(flens_INCLUDE_DIRS ${install_dir}/include)
  set(flens_FOUND 1)
endif()

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
#! ### Installation of dependencies
#! <!--#########################################################-->
#!
#! If blitz and flens are compiled as external projects, install them as well. Libraries
#! and header files are installed into the same locations as CPPQEDcore, this way it is
#! ensured they are found.

set(BUNDLED_HEADERS "${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_DIR}/bundled")

if(BUNDLED_BLITZ)
  ExternalProject_Get_Property(cppqed-blitz install_dir)
  file(GLOB blitz_libs ${install_dir}/lib/*.so*)
  set(blitz_libs ${install_dir}/lib/libblitzcppqed.so
                 ${install_dir}/lib/libblitzcppqed.so.0
                 ${install_dir}/lib/libblitzcppqed.so.0.0.0)
  install(FILES ${blitz_libs}
          PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${install_dir}/include/blitz ${install_dir}/include/random
          DESTINATION ${BUNDLED_HEADERS})
endif()

if(BUNDLED_FLENS)
  ExternalProject_Get_Property(flens install_dir)
  install(DIRECTORY ${install_dir}/include/flens
          DESTINATION ${BUNDLED_HEADERS})
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
