# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

cmake_minimum_required (VERSION 3.12)

# ### Version management
#
# Set the major, minor and patch version of C++QED (c.f. \ref versioning). The major and minor
# versions end up in the library names of all components. In this section,
# also the ABI version of the library is set (c.f. \ref ABI).

# MAJOR - the grand version of C++QED, as in v2
set(CPPQED_VERSION_MAJOR 3)
# MINOR - the milestone release, as in v2m10
set(CPPQED_VERSION_MINOR 0)
# PATCH - should be increased on every snapshot package release, reset on milestone change
set(CPPQED_VERSION_PATCH 0)
# Full C++QED version
set(CPPQED_VERSION
  "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}.${CPPQED_VERSION_PATCH}")
# MAJOR.MINOR, at the moment mainly distinguishes stable (2.10) from development (2.100)
set(CPPQED_ID "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}")

project(cppqed VERSION ${CPPQED_VERSION}  
               LANGUAGES CXX
               DESCRIPTION "A framework for simulating open quantum dynamics"
               HOMEPAGE_URL "https://github.com/vukics/cppqed")

include(CTest)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# otherwise object libraries are not compiled with the -fPIC option
# cf. http://gitlab.kitware.com/cmake/cmake/-/issues/17356

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

option(ENABLE_TESTING "Enable testsuite" ON)

if(ENABLE_TESTING)
  enable_testing()
endif()

# Build type: guard against values other than debug/release, make release default

if (NOT CMAKE_BUILD_TYPE)
  message(WARNING "Build type not set, default is \"Release\".")
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type_tolower)
if (NOT(cmake_build_type_tolower MATCHES "^(release|debug)$"))
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug and Release (case-insensitive).")
endif()

set(CMAKE_DEBUG_POSTFIX "_d")

# ### Compiler detection
# Ubuntu 20.04LTS versions

set(GCC_MINIMAL 9.3)
set(CLANG_MINIMAL 10.0)
set(BOOST_MINIMAL 1.71)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${GCC_MINIMAL})
  message(FATAL_ERROR "GCC g++ version >= ${GCC_MINIMAL} needed.")
endif ()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MINIMAL})
  message(FATAL_ERROR "Clang version >= ${CLANG_MINIMAL} needed.")
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

include(GetGitRevisionDescription)
get_git_head_revision(CPPQED_GIT_REFSPEC CPPQED_GIT_SHA1)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppqedGitSHA1.cc.in" "${CMAKE_CURRENT_BINARY_DIR}/cppqedGitSHA1.cc" @ONLY)

add_library(cppqedGitSHA1 OBJECT "${CMAKE_CURRENT_BINARY_DIR}/cppqedGitSHA1.cc")

add_compile_options(-Wall -Wextra -Wpedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-ignored-qualifiers 
-Wno-sign-compare -Wno-overloaded-virtual -Wunreachable-code -Wno-unknown-pragmas -Wno-narrowing)

#macro(add_executable S)
  #if (${cmake_build_type_tolower} STREQUAL debug)
    #_add_executable("${S}_d" ${ARGN})
  #else()
    #_add_executable(${S} ${ARGN})
  #endif()
#endmacro()

add_subdirectory(CPPQEDutils)

add_subdirectory(CPPQEDcore)

add_subdirectory(CPPQEDelements)

add_subdirectory(CPPQEDscripts)

add_subdirectory(cpypyqed)

if(ENABLE_TESTING)
  add_subdirectory(Testing)
endif()

if (EXISTS ${PROJECT_SOURCE_DIR}/OwnScripts)
  add_subdirectory(OwnScripts)
endif()
