# - Try to find blitz++
# Once done, this will define
#
#  blitz++_FOUND - system has blitz++
#  blitz++_INCLUDE_DIRS - the blitz++ include directories
#  blitz++_LIBRARIES - link these to use blitz++

include(LibFindMacros)
include (CheckIncludeFiles) 

# Include dir
find_path(blitz++_INCLUDE_DIR
  NAMES blitz/tinyvec2.h
)

# Finally the library itself
find_library(blitz++_LIBRARY
  NAMES blitz
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(blitz++_PROCESS_INCLUDES blitz++_INCLUDE_DIR)
set(blitz++_PROCESS_LIBS blitz++_LIBRARY)
libfind_process(blitz++)

file( STRINGS ${blitz++_INCLUDE_DIR}/blitz/gnu/bzconfig.h BLITZ_SERIALIZATION REGEX "^#define BZ_HAVE_BOOST_SERIALIZATION")