# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

cmake_minimum_required (VERSION 3.12)

# ### Version management
#
# Set the major, minor and patch version of C++QED (c.f. \ref versioning). The major and minor
# versions end up in the library names of all components. In this section,
# also the ABI version of the library is set (c.f. \ref ABI).

# MAJOR - the grand version of C++QED, as in v2
set(CPPQED_VERSION_MAJOR 2)
# MINOR - the milestone release, as in v2m10
set(CPPQED_VERSION_MINOR 120)
# PATCH - should be increased on every snapshot package release, reset on milestone change
set(CPPQED_VERSION_PATCH 0)
# Full C++QED version
set(CPPQED_VERSION
  "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}.${CPPQED_VERSION_PATCH}")
# MAJOR.MINOR, at the moment mainly distinguishes stable (2.10) from development (2.100)
set(CPPQED_ID "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}")

# ABI versioning scheme
#
# Adopt the shared library versioning scheme of libtool
# - CURRENT  - the current ABI version
# - AGE      - number of versions backwards that CURRENT is compatible with
# - REVISION - patch revision
#
# *Rationale*
#
# Bugfixes and Patches which do not affect ABI:
#  * increase REVISION
#
# Symbols added to library (i.e. binary compatibility NOT broken):
#  * increase AGE
#  * remember to call dh_makeshlibs with -V packagename (>=packageversion) when packaging
#
# Symbols removed or changed (i.e. binary compatibility broken):
#  * set CURRENT -> CURRENT + AGE + 1, reset AGE
#
# On every change of AGE or CURRENT: reset REVISION
#
# The library version (.so.X.Y.Z) is `{CURRENT-AGE}.{AGE}.{REVISION}`, here named
# `{ABI_MAJOR}.{ABI_MINOR}.{ABI_MICRO}`. This way it is ensured that ABI_MAJOR only changes
# when binary compatibility is broken.
#
# The SONAME of the library always is: `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}.so.${ABI_MAJOR}`
# and the packages are named `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}-${ABI_MAJOR}[-dev]`

# The current ABI version
set(CPPQED_ABI_CURRENT 6)
# number of Versions backwards that CURRENT is compatible with
set(CPPQED_ABI_AGE 0)
# Patch revision
set(CPPQED_ABI_REVISION 0)
# Set to ABI_CURRENT - ABI_AGE
set(CPPQED_ABI_MAJOR)
math(EXPR CPPQED_ABI_MAJOR "${CPPQED_ABI_CURRENT}-${CPPQED_ABI_AGE}")
# Same as ABI_AGE
set(CPPQED_ABI_MINOR "${CPPQED_ABI_AGE}")
# Same as ABI_REVISION
set(CPPQED_ABI_MICRO "${CPPQED_ABI_REVISION}")

project(cppqed VERSION ${CPPQED_VERSION}  
               LANGUAGES CXX
               DESCRIPTION "A framework for simulating open quantum dynamics"
               HOMEPAGE_URL "https://github.com/vukics/cppqed")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

# Build type: guard against values other than debug/release, make release default

if (NOT CMAKE_BUILD_TYPE)
  message(WARNING "Build type not set, default is \"Release\".")
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type_tolower)
if(NOT(cmake_build_type_tolower MATCHES "^(release|debug)$"))
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug and Release (case-insensitive).")
endif()

set(CMAKE_DEBUG_POSTFIX "_d")

# ### Compiler detection
# Ubuntu 20.04LTS versions

set(GCC_MINIMAL 9.3)
set(CLANG_MINIMAL 10.0)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${GCC_MINIMAL})
  message(FATAL_ERROR "GCC g++ version >= ${GCC_MINIMAL} needed.")
endif ()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MINIMAL})
  message(FATAL_ERROR "Clang version >= ${CLANG_MINIMAL} needed.")
endif ()

# Switch for FLENS support.
# option(FLENS "FLENS support" ON)

message(STATUS "Configuring blitz")
add_subdirectory(blitz)
include_directories(SYSTEM blitz) # in order that we do not receive warnings from blitz

add_compile_options(-Wall -Wextra -Wpedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-ignored-qualifiers -Wno-sign-compare -Wno-overloaded-virtual)

add_subdirectory(CPPQEDutils)

# add_subdirectory(CPPQEDcore)

