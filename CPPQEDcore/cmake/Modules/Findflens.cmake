# Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \file
#! \ingroup FindPackage
#! \brief Find FLENS.
#!
#! Once done, this will define
#!
#! - `flens_FOUND` - system has flens
#! - `flens_INCLUDE_DIRS` - the flens include directories

if (UNIX)
  find_package(PkgConfig QUIET)
  pkg_check_modules(flens_PKGCONF QUIET flens)
endif()

# Include dir
find_path(flens_INCLUDE_DIR
  NAMES flens/flens.cxx
  PATHS ${flens_PKGCONF_INCLUDE_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(flens_PROCESS_INCLUDES flens_INCLUDE_DIR)
libfind_process(flens)