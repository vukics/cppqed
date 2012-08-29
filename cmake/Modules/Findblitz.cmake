# - Try to find blitz
# Once done, this will define
#
#  blitz_FOUND - system has blitz
#  blitz_INCLUDE_DIRS - the blitz include directories
#  blitz_LIBRARIES - link these to use blitz
#  blitz_SERIALIZATION_FOUND - true if blitz was configured with --enable-serialization

include(LibFindMacros)
include(CheckIncludeFiles)

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

if(${CMAKE_VERSION} VERSION_GREATER "2.8.5")
  # the symbol blitz::paddedArray is new in blitz++-0.10,we check if it exists in the library
  # if not we might have accidentially picked up an older version of the library.
  include(CheckCXXSymbolExists)
  set(CMAKE_REQUIRED_LIBRARIES ${blitz_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${blitz_INCLUDE_DIRS})
  CHECK_CXX_SYMBOL_EXISTS(blitz::paddedArray ${blitz_INCLUDE_DIR}/blitz/array.h blitz_IS_MERCURIAL_VERSION)
  set(CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_INCLUDES)
  if(NOT blitz_IS_MERCURIAL_VERSION)
    set(blitz_LIBRARY)
    message(FATAL_ERROR "Blitz++ version >= 0.10 (or mercurial checkout) is needed.")
  endif(NOT blitz_IS_MERCURIAL_VERSION)
  CHECK_SYMBOL_EXISTS(BZ_HAVE_BOOST_SERIALIZATION ${blitz_INCLUDE_DIR}/blitz/gnu/bzconfig.h blitz_SERIALIZATION_FOUND)
else(${CMAKE_VERSION} VERSION_GREATER "2.8.5")
  file( STRINGS ${blitz_INCLUDE_DIR}/blitz/gnu/bzconfig.h blitz_SERIALIZATION_FOUND REGEX "^#define BZ_HAVE_BOOST_SERIALIZATION")
endif(${CMAKE_VERSION} VERSION_GREATER "2.8.5")

