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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#! \file
#! <!--#########################################################-->
#! ### Version management
#! <!--#########################################################-->
#!
#! Set the major, minor and patch version of C++QED (c.f. \ref versioning). The major and minor
#! versions end up in the library names of all components. In this section,
#! also the ABI version of the library is set (c.f. \ref ABI).


#! \anchor versioning
#! \name Version variables
#! Versioning scheme
#! @{

#! MAJOR - the grand version of C++QED, as in v2
set(CPPQED_VERSION_MAJOR 2)
#! MINOR - the milestone release, as in v2m10
set(CPPQED_VERSION_MINOR 120)
#! PATCH - should be increased on every snapshot package release, reset on milestone change
set(CPPQED_VERSION_PATCH 10)
#! Full C++QED version
set(CPPQED_VERSION
  "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}.${CPPQED_VERSION_PATCH}")
#! MAJOR.MINOR, at the moment mainly distinguishes stable (2.10) from development (2.100)
set(CPPQED_ID "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}")

#! @}

#! \anchor ABI
#! \name ABI versioning scheme
#!
#! Adopt the shared library versioning scheme of libtool
#! - CURRENT  - the current ABI version
#! - AGE      - number of versions backwards that CURRENT is compatible with
#! - REVISION - patch revision
#!
#! *Rationale*
#!
#! Bugfixes and Patches which do not affect ABI:
#!  * increase REVISION
#!
#! Symbols added to library (i.e. binary compatibility NOT broken):
#!  * increase AGE
#!  * remember to call dh_makeshlibs with -V packagename (>=packageversion) when packaging
#!
#! Symbols removed or changed (i.e. binary compatibility broken):
#!  * set CURRENT -> CURRENT + AGE + 1, reset AGE
#!
#! On every change of AGE or CURRENT: reset REVISION
#!
#! The library version (.so.X.Y.Z) is `{CURRENT-AGE}.{AGE}.{REVISION}`, here named
#! `{ABI_MAJOR}.{ABI_MINOR}.{ABI_MICRO}`. This way it is ensured that ABI_MAJOR only changes
#! when binary compatibility is broken.
#!
#! The SONAME of the library always is: `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}.so.${ABI_MAJOR}`
#! and the packages are named `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}-${ABI_MAJOR}[-dev]`
#! @{

#! The current ABI version
set(CPPQED_ABI_CURRENT 5)
#! number of Versions backwards that CURRENT is compatible with
set(CPPQED_ABI_AGE 0)
#! Patch revision
set(CPPQED_ABI_REVISION 0)
#! Set to ABI_CURRENT - ABI_AGE
set(CPPQED_ABI_MAJOR)
math(EXPR CPPQED_ABI_MAJOR "${CPPQED_ABI_CURRENT}-${CPPQED_ABI_AGE}")
#! Same as ABI_AGE
set(CPPQED_ABI_MINOR "${CPPQED_ABI_AGE}")
#! Same as ABI_REVISION
set(CPPQED_ABI_MICRO "${CPPQED_ABI_REVISION}")

#! @}

#! \file
#! <!--#########################################################-->
#! ### Compiler detection
#! <!--#########################################################-->
#!
#! Ubuntu 20.04LTS versions

set(GCC_MINIMAL 9.3)
set(CLANG_MINIMAL 10.0)


if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${GCC_MINIMAL})
  message(FATAL_ERROR "GCC g++ version >= ${GCC_MINIMAL} needed.")
endif ()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MINIMAL})
  message(FATAL_ERROR "Clang version >= ${CLANG_MINIMAL} needed.")
endif ()

option(BUNDLED_BLITZ "Download and compile blitz++ automatically" OFF)

#! Switch for FLENS support.
option(FLENS "FLENS support" ON)

#! If this option is set, %CMake downloads and compiles the optional FLENS dependency automatically.
option(BUNDLED_FLENS "Download and compile FLENS automatically" OFF)

if(BUNDLED_BLITZ OR BUNDLED_FLENS)
  find_package(Git 1.6.5 REQUIRED) # ExternalProject requires this version
endif()


if(BUNDLED_BLITZ)
  include(ExternalProject)

#   if( CMAKE_SIZEOF_VOID_P EQUAL 8 ) # 64bit system
#     set(BLITZ_CONFIG_OPTIONS "${BLITZ_CONFIG_OPTIONS} --enable-64bit ")
#   else( CMAKE_SIZEOF_VOID_P EQUAL 8 ) # 32bit system
#     set(BLITZ_CONFIG_OPTIONS "${BLITZ_CONFIG_OPTIONS} --enable-simd-width=8 ")
#   endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

  if(EXISTS ${CMAKE_SOURCE_DIR}/external/blitz) # distributed packaged together with cppqed
    set(BLITZ_SOURCE DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/blitz/. <SOURCE_DIR>/)
  else()
    set(BLITZ_SOURCE GIT_REPOSITORY http://github.com/vukics/blitz.git GIT_TAG a39b77e33a92cce48495306817dc1573f62f755e GIT_SHALLOW TRUE GIT_PROGRESS TRUE)
  endif()

  ExternalProject_Add(
    blitz
    ${BLITZ_SOURCE}
    UPDATE_COMMAND ""
    CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  )
  ExternalProject_Get_Property(blitz install_dir)
  message(STATUS "Bundled blitz to be installed to ${install_dir}")
  
  set(blitz_INCLUDE_DIRS ${install_dir}/include)
  set(blitz_LIBRARIES ${install_dir}/lib/libblitz${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(blitz_SERIALIZATION_FOUND 1)
  set(blitz_FOUND 1)
endif()

add_subdirectory(CPPQEDutils)

