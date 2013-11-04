cmake_minimum_required (VERSION 2.8.9)

project (CPPQEDcore)

get_filename_component(CPPQED_CMAKE_MODULE_PATH "cmake/Modules" ABSOLUTE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CPPQED_CMAKE_MODULE_PATH})

include(FeatureSummary)
include(CPPQEDUse)

##################################################
# Version
##################################################

# Versioning scheme

# MAJOR - the grand version of C++QED, as in v2
# MINOR - the milestone release, as in v2m9
# PATCH - should be increased on every snapshot package release, reset on milestone change

set(CPPQED_MAJOR_VERSION 2)
set(CPPQED_MINOR_VERSION 99)
set(CPPQED_PATCH_VERSION 0)
set(CPPQED_VERSION
  ${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION}.${CPPQED_PATCH_VERSION})
set(CPPQED_ID ${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION})


# Adopt the shared library versioning scheme of libtool
# CURRENT  - the current ABI version
# AGE      - number of versions backwards that CURRENT is compatible with
# REVISION - patch revision

# Rationale:
# Bugfixes and Patches which do not affect ABI:
#  * increase REVISION
# Symbols added to library (i.e. binary compatibility NOT broken):
#  * increase AGE
#  * remember to call dh_makeshlibs with -V packagename (>=packageversion) when packaging
# Symbols removed or changed (i.e. binary compatibility broken):
#  * set CURRENT -> CURRENT + AGE + 1, reset AGE
# On every change of AGE or CURRENT: reset REVISION

# The library version (.so.X.Y.Z) is {CURRENT-AGE}.{AGE}.{REVISION}, here named
# {ABI_MAJOR}.{ABI_MINOR}.{ABI_MICRO}. This way it is ensured that ABI_MAJOR only changes
# when binary compatibility is broken.

# The SONAME of the library always is: libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}.so.${ABI_MAJOR}
# and the packages are named libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}-${ABI_MAJOR}[-dev]

set(CPPQED_ABI_CURRENT 0)
set(CPPQED_ABI_AGE 0)
set(CPPQED_ABI_REVISION 0)
math(EXPR CPPQED_ABI_MAJOR ${CPPQED_ABI_CURRENT}-${CPPQED_ABI_AGE})
set(CPPQED_ABI_MINOR ${CPPQED_ABI_AGE})
set(CPPQED_ABI_MICRO ${CPPQED_ABI_REVISION})
  
##################################################
# Installation directories
##################################################

set(CPPQED_CMAKE_SUBDIR "cmake/CPPQED-${CPPQED_ID}")
set(CPPQED_INCLUDE_SUBDIR "CPPQED-${CPPQED_ID}/core")

##################################################
# Compiler detection
##################################################
set(G++_MINIMAL 4.7)
set(CLANG_MINIMAL 3.1)


if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${G++_MINIMAL})
  message(FATAL_ERROR "GCC g++ version >= ${G++_MINIMAL} needed.")
endif ()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MINIMAL})
  message(FATAL_ERROR "Clang version >= ${CLANG_MINIMAL} needed.")
endif ()


##################################################
# Library detection
##################################################

find_package(GSL REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

# The option SERIALIZATION will be disabled if any of the compoments required for serialization is missing
option(SERIALIZATION "Boost serialization" ON)
if(SERIALIZATION)
  set(CPPQED_HAS_SERIALIZATION 1)
else(SERIALIZATION)
  set(CPPQED_HAS_SERIALIZATION 0)
endif(SERIALIZATION)

# The option FLENS will be disabled if flens is missing
option(FLENS "FLENS support" ON)

# blitz
find_package(blitz REQUIRED)
if( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built with serialization support.")
else( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built without serialization support. Please configure blitz with --enable-serialization to enable optional binary statevector output.")
  set(CPPQED_HAS_SERIALIZATION 0)
endif( blitz_SERIALIZATION_FOUND )

include_directories(SYSTEM ${blitz_INCLUDE_DIRS} )
set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS})

# Boost
if(${CMAKE_VERSION} VERSION_LESS "2.8.8")
  find_package(Boost REQUIRED)
  find_package(Boost COMPONENTS serialization)
else(${CMAKE_VERSION} VERSION_LESS "2.8.8")
  find_package(Boost REQUIRED OPTIONAL_COMPONENTS serialization)
endif(${CMAKE_VERSION} VERSION_LESS "2.8.8")

if( Boost_SERIALIZATION_FOUND )
  message(STATUS "Boost serialization library found.")
else( Boost_SERIALIZATION_FOUND )
  message(STATUS "Boost serialization library not found.")
  set(CPPQED_HAS_SERIALIZATION 0)
endif( Boost_SERIALIZATION_FOUND )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} )
set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})


# FLENS
find_package(flens)
if( flens_FOUND AND FLENS )
  include_directories(SYSTEM ${flens_INCLUDE_DIRS})
  #set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${flens_INCLUDE_DIRS})
  set(PUBLIC_LIBS ${PUBLIC_LIBS} ${flens_LIBRARIES})
  set(CPPQED_HAS_FLENS 1)
else( flens_FOUND AND FLENS )
  message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  set(DO_NOT_USE_FLENS ON)
  set(CPPQED_HAS_FLENS 0)
endif( flens_FOUND AND FLENS )

# Check if serialization can be enabled, inform user
if( CPPQED_HAS_SERIALIZATION )
  message(STATUS "Support for binary statevector output enabled.")
  set(PUBLIC_LIBS ${PUBLIC_LIBS} ${Boost_SERIALIZATION_LIBRARY})
else( CPPQED_HAS_SERIALIZATION )
  set(DO_NOT_USE_BOOST_SERIALIZATION ON)
  message(STATUS "Optional support for binary statevector output disabled.")
endif( CPPQED_HAS_SERIALIZATION )

set(PRIVATE_LIBS ${PRIVATE_LIBS} ${GSL_LIBRARIES})
set(PUBLIC_LIBS ${PUBLIC_LIBS} ${blitz_LIBRARIES})

##################################################
# Flags for release and debug mode
##################################################

# Generate config.h
configure_file(cppqedconfig.h.in cppqedconfig.h)

if(Boost_FOUND)
  if(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    set(CPPQED_DEFINITIONS ${CPPQED_DEFINITIONS} -DBOOST_RESULT_OF_USE_TR1)
    message(STATUS "added -DBOOST_RESULT_OF_USE_TR1" )
  endif()
endif()

if(CBLAS_FOUND AND "${CBLAS_LIBRARIES}" MATCHES "gslcblas")
  set(CPPQED_DEFINITIONS ${CPPQED_DEFINITIONS} -DGSL_CBLAS)
  message(STATUS "added -DGSL_CBLAS" )
endif(CBLAS_FOUND AND "${CBLAS_LIBRARIES}" MATCHES "gslcblas")

##################################################
# Compilation
##################################################
CPPQED_SETUP()

include_directories(${PROJECT_BINARY_DIR}) # for cppqedconfig.h

# build libC++QEDcore.so
set(CORE_SOURCE_DIRS utils quantumdata quantumoperator structure quantumtrajectory composites)
set(quantumdata_NEEDS utils)
set(quantumoperator_NEEDS quantumdata utils)
set(structure_NEEDS quantumdata utils)
set(quantumtrajectory_NEEDS structure quantumdata utils)
set(composites_NEEDS structure quantumdata utils)
foreach(d ${CORE_SOURCE_DIRS})
  add_subdirectory(${d})
endforeach(d)
set(CPPQEDLIB C++QED-${CPPQED_ID})

gather_includes(CPPQED CORE_SOURCE_DIRS)
add_library(${CPPQEDLIB} SHARED ${OBJ_TARGETS} ${CPPQED_PUBLIC_HEADERS} cppqedconfig.h.in)
set(CPPQED_PUBLIC_HEADERS ${CPPQED_PUBLIC_HEADERS} ${PROJECT_BINARY_DIR}/cppqedconfig.h)
target_link_libraries(${CPPQEDLIB} LINK_PUBLIC ${PUBLIC_LIBS} LINK_PRIVATE ${PRIVATE_LIBS})
set_target_properties(${CPPQEDLIB} PROPERTIES
      PUBLIC_HEADER "${CPPQED_PUBLIC_HEADERS}"
      INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
      VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
      SOVERSION ${CPPQED_ABI_MAJOR}
)
##################################################
# Installation
##################################################

install(TARGETS ${CPPQEDLIB}
        EXPORT CPPQEDcoreTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR}
        COMPONENT shlib
)

# Add all targets to the build-tree export set
export(TARGETS ${CPPQEDLIB}
  FILE "${PROJECT_BINARY_DIR}/CPPQEDcoreTargets.cmake")
export(PACKAGE CPPQED)
  
# Create the CPPQEDConfig.cmake
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_BINARY_DIR}")
foreach(d ${CORE_SOURCE_DIRS})
  set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/${d}) 
endforeach(d)
set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
set(CONF_FLAT_INCLUDE 0)
configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}/CPPQEDConfig.cmake"
  INSTALL_DESTINATION  "${PROJECT_BINARY_DIR}"
  PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR
)
write_basic_package_version_file(${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake 
  VERSION ${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION}
  COMPATIBILITY ExactVersion
)
configure_file(cmake/Modules/CPPQEDUse.cmake ${PROJECT_BINARY_DIR}/CPPQEDUse.cmake COPYONLY)
configure_file(cmake/Modules/ElementsTemplateConfig.cmake.in ${PROJECT_BINARY_DIR}/ElementsTemplateConfig.cmake.in COPYONLY)

# ... and for the installation tree
set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR})
set(CONF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR})
set(CONF_FLAT_INCLUDE 1)
configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  INSTALL_DESTINATION "${CONF_CMAKE_DIR}"
  PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR
)

# Install the CPPQEDConfig.cmake and CPPQEDConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  "${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}/CPPQEDUse.cmake"
  "${PROJECT_BINARY_DIR}/ElementsTemplateConfig.cmake.in"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT CPPQEDcoreTargets DESTINATION
  "${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)

  
##################################################
# Summary of enabled/disabled features
##################################################

set_package_properties(PkgConfig PROPERTIES URL "http://pkgconfig.freedesktop.org/wiki"
                                DESCRIPTION "Package config system that manages compile/link flags"
                                TYPE OPTIONAL
                                PURPOSE "Assists cmake in finding libraries.")
set_package_properties(blitz PROPERTIES  URL "http://sf.net/projects/blitz/"
                                DESCRIPTION "High-performance C++ vector mathematics library"
                                TYPE REQUIRED
                                PURPOSE "Multi-Array implementation used in the framework.")
set_package_properties(flens PROPERTIES URL "http://www.mathematik.uni-ulm.de/~lehn/FLENS/"
                                DESCRIPTION "Flexible Library for Efficient Numerical Solutions."
                                TYPE OPTIONAL
                                PURPOSE "With FLENS the partial transpose of the density operator can be calculated.")
set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                DESCRIPTION "Collection of portable C++ source libraries."
                                TYPE REQUIRED
                                PURPOSE "Advanced template metaprogramming and preprocessor algorithms.")
set_package_properties(GSL PROPERTIES URL "http://www.gnu.org/software/gsl/"
                                DESCRIPTION "GNU Scientific Library"
                                TYPE REQUIRED
                                PURPOSE "Used as implementation of an ODE solver and random number generator.")
set_package_properties(CBLAS PROPERTIES URL "http://math-atlas.sourceforge.net/"
                                DESCRIPTION "C library for Basic Linear Algebra Subprograms"
                                TYPE OPTIONAL
                                PURPOSE "Needed for GSL.")
add_feature_info("FLENS" CPPQED_HAS_FLENS "compile framework with FLENS support.")
add_feature_info(Serialization CPPQED_HAS_SERIALIZATION "needed for binary statevector output." )
feature_summary( WHAT ALL )
add_subdirectory(examples)