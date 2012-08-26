# - Find CBLAS
# Find the native CBLAS headers and libraries.
#
#  CBLAS_LIBRARIES    - List of libraries when using cblas.
#  CBLAS_INCLUDE_DIRS - List of include directories
#  CBLAS_FOUND        - True if cblas found.
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
include(LibFindMacros)

# Include dir
set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
find_library(CBLAS_LIBRARY
  NAMES cblas gslcblas
  PATHS $ENV{CBLASDIR}/lib $ENV{CBLASDIR}/lib64 $ENV{UIBK_GSL_LIB}
)

if(${CBLAS_LIBRARY} MATCHES gslcblas)
  set(CBLAS_INCLUDE_CANDIDATE gsl/gsl_cblas.h)
else(${CBLAS_LIBRARY} MATCHES gslcblas)
  set(CBLAS_INCLUDE_CANDIDATE cblas.h)
endif(${CBLAS_LIBRARY} MATCHES gslcblas)

find_path(CBLAS_INCLUDE_DIR ${CBLAS_INCLUDE_CANDIDATE} HINTS $ENV{CBLASDIR}/include $ENV{UIBK_GSL_INC})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CBLAS_PROCESS_INCLUDES CBLAS_INCLUDE_DIR)
set(CBLAS_PROCESS_LIBS CBLAS_LIBRARY)
libfind_process(CBLAS)
