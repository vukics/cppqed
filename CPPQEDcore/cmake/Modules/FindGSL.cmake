# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

#! \file
#! \ingroup FindPackage
#! \brief Find GSL.
#!
#! Find the native GSL headers and libraries.
#!
#! - `GSL_INCLUDE_DIRS` - where to find gsl/gsl_linalg.h, etc.
#! - `GSL_LIBRARIES` - List of libraries when using gsl.
#! - `GSL_FOUND` - True if gsl found.
#!

include(LibFindMacros)
libfind_package(GSL CBLAS)

if (UNIX)
  find_package(PkgConfig QUIET)
  pkg_check_modules(GSL_PKGCONF QUIET gsl)
endif()

find_path(GSL_INCLUDE_DIR gsl/gsl_linalg.h HINTS ${GSL_PKGCONF_INCLUDE_DIRS} $ENV{UIBK_GSL_INC})
find_library(GSL_LIBRARY NAMES gsl HINTS ${GSL_PKGCONF_LIBRARY_DIRS} $ENV{UIBK_GSL_LIB})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(GSL_PROCESS_INCLUDES GSL_INCLUDE_DIR CBLAS_INCLUDE_DIRS)
set(GSL_PROCESS_LIBS GSL_LIBRARY CBLAS_LIBRARIES)
libfind_process(GSL)
