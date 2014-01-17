#! \file
#! \ingroup FindPackage
#! \brief Find FLENS.
#!
#! Once done, this will define
#!
#! - `flens_FOUND` - system has flens
#! - `flens_INCLUDE_DIRS` - the flens include directories
#! - `flens_LIBRARIES` - link these to use flens

include(LibFindMacros)
libfind_package(flens CBLAS)

if (UNIX)
  find_package(PkgConfig QUIET)
  pkg_check_modules(flens_PKGCONF QUIET flens)
endif()

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