# Copyright Raimar Sandner 2012–2020. Copyright András Vukics 2020-2021. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \file
#! \ingroup FindPackage
#! \brief Find blitz++.
#!
#! Once done, this will define
#!
#! - `blitz_FOUND`: system has blitz
#! - `blitz_INCLUDE_DIRS`: the blitz include directories
#! - `blitz_LIBRARIES`: link these to use blitz
#! - `blitz_SERIALIZATION_FOUND`: true if blitz was configured with `--enable-serialization`
#!
#! Note that blitz will only be considered found if it is the [C++QED patched version](https://github.com/vukics/blitz.git).
#! This is checked by looking if it has `BLITZ_ARRAY_LARGEST_RANK`.

include(LibFindMacros)
include(CheckCXXSymbolExists)

if (UNIX)
  find_package(PkgConfig QUIET)
  pkg_check_modules(blitz_PKGCONF QUIET blitz)
endif()

# Include dir
find_path(blitz_INCLUDE_DIR
  NAMES blitz/tinyvec2.h
  HINTS ${blitz_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(blitz_LIBRARY
  NAMES blitz
  HINTS ${blitz_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(blitz_PROCESS_INCLUDES blitz_INCLUDE_DIR)
set(blitz_PROCESS_LIBS blitz_LIBRARY)
libfind_process(blitz)

CHECK_CXX_SYMBOL_EXISTS(BLITZ_ARRAY_LARGEST_RANK "${blitz_INCLUDE_DIR}/blitz/blitz.h" blitz_IS_CPPQED_VERSION)
if(NOT blitz_IS_CPPQED_VERSION)
  message(FATAL_ERROR "You need the C++QED version of blitz++. The repository is at https://github.com/vukics/blitz.git")
endif(NOT blitz_IS_CPPQED_VERSION)

file(GLOB BZCONFIG ${blitz_INCLUDE_DIR}/blitz/*/bzconfig.h)
CHECK_CXX_SYMBOL_EXISTS(BZ_HAVE_BOOST_SERIALIZATION "${BZCONFIG}" blitz_SERIALIZATION_FOUND)

if (blitz_FOUND)
  if (NOT blitz_FIND_QUIETLY)
    message (STATUS "Found components for blitz++: includes folder = ${blitz_INCLUDE_DIRS} ; libraries = ${blitz_LIBRARIES}")
  endif (NOT blitz_FIND_QUIETLY)
else (blitz_FOUND)
  if (blitz_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find Blitz++!")
  endif (blitz_FIND_REQUIRED)
endif (blitz_FOUND)
