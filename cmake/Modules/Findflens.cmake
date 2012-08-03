# - Try to find flens
# Once done, this will define
#
#  flens_FOUND - system has flens
#  flens_INCLUDE_DIRS - the flens include directories
#  flens_LIBRARIES - link these to use flens

include(LibFindMacros)

# Include dir
find_path(flens_INCLUDE_DIR
  NAMES flens/flens.h
)

# Finally the library itself
find_library(blitz++_LIBRARY
  NAMES flens
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(flens_PROCESS_INCLUDES flens_INCLUDE_DIR)
set(flens_PROCESS_LIBS flens_LIBRARY)
libfind_process(flens)