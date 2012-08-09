# - Try to find flens
# Once done, this will define
#
#  flens_FOUND - system has flens
#  flens_INCLUDE_DIRS - the flens include directories
#  flens_LIBRARIES - link these to use flens
#  flens_DEFINITIONS - adds -DGSL_CBLAS if gslcblas is detected

include(LibFindMacros)
libfind_package(flens CBLAS)
if(CBLAS_FOUND AND ${CBLAS_LIBRARY} MATCHES "gslcblas")
  set(flens_DEFINITIONS "-DGSL_CBLAS")
else(CBLAS_FOUND AND ${CBLAS_LIBRARY} MATCHES "gslcblas")
  set(flens_DEFINITIONS "")
endif(CBLAS_FOUND AND ${CBLAS_LIBRARY} MATCHES "gslcblas")

libfind_pkg_check_modules(flens_PKGCONF flens)

# Include dir
find_path(flens_INCLUDE_DIR
  NAMES flens/flens.h
  PATHS ${flens_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(flens_LIBRARY
  NAMES flens
  PATHS ${flens_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(flens_PROCESS_INCLUDES flens_INCLUDE_DIR CBLAS_INCLUDE_DIRS)
set(flens_PROCESS_LIBS flens_LIBRARY CBLAS_LIBRARIES)
libfind_process(flens)