# - Try to find blitz++
# Once done, this will define
#
#  blitz++_FOUND - system has blitz++
#  blitz++_INCLUDE_DIRS - the blitz++ include directories
#  blitz++_LIBRARIES - link these to use blitz++
#  blitz++_SERIALIZATION_FOUND - true if blitz++ was configured with --enable-serialization

include(LibFindMacros)
include (CheckIncludeFiles) 

libfind_pkg_check_modules(blitz++_PKGCONF blitz++)

# Include dir
find_path(blitz++_INCLUDE_DIR
  NAMES blitz/tinyvec2.h
  PATHS ${blitz++_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(blitz++_LIBRARY
  NAMES blitz
  PATHS ${blitz++_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(blitz++_PROCESS_INCLUDES blitz++_INCLUDE_DIR)
set(blitz++_PROCESS_LIBS blitz++_LIBRARY)
libfind_process(blitz++)

file( STRINGS ${blitz++_INCLUDE_DIR}/blitz/gnu/bzconfig.h BLITZ_SERIALIZATION REGEX "^#define BZ_HAVE_BOOST_SERIALIZATION")

if( BLITZ_SERIALIZATION )
  set(blitz++_SERIALIZATION_FOUND 1)
else( BLITZ_SERIALIZATION )
  set(blitz++_SERIALIZATION_FOUND 0)
endif( BLITZ_SERIALIZATION)
libfind_process(blitz++_SERIALIZATION)