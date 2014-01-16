
#! \file
#! \ingroup FindPackage
#! \brief Try to find blitz.
#!
#! Once done, this will define
#!
#! - `blitz_FOUND`: system has blitz
#! - `blitz_INCLUDE_DIRS`: the blitz include directories
#! - `blitz_LIBRARIES`: link these to use blitz
#! - `blitz_SERIALIZATION_FOUND`: true if blitz was configured with `--enable-serialization`
#!
#! Not that blitz will only be considered found if it is the [C++QED patched version](http://sourceforge.net/p/cppqed/blitz/ci/default/tree/).
#! This is checked by looking if it has `BLITZ_ARRAY_LARGEST_RANK`.

include(LibFindMacros)
include(CheckIncludeFiles)
include(CheckSymbolExists)

libfind_pkg_check_modules(blitz_PKGCONF blitz)

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
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(blitz_PROCESS_INCLUDES blitz_INCLUDE_DIR)
set(blitz_PROCESS_LIBS blitz_LIBRARY)
libfind_process(blitz)

file(READ ${blitz_INCLUDE_DIR}/blitz/blitz.h BLITZ_H)
string(REGEX MATCH BLITZ_ARRAY_LARGEST_RANK blitz_IS_CPPQED_VERSION ${BLITZ_H})
if(NOT blitz_IS_CPPQED_VERSION)
  message(FATAL_ERROR "You need the C++QED version of blitz++. The repository is at http://sourceforge.net/p/cppqed/blitz/ci/default/tree/")
endif(NOT blitz_IS_CPPQED_VERSION)

file(GLOB BZCONFIG ${blitz_INCLUDE_DIR}/blitz/*/bzconfig.h)
CHECK_SYMBOL_EXISTS(BZ_HAVE_BOOST_SERIALIZATION ${BZCONFIG} blitz_SERIALIZATION_FOUND)
