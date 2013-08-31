cmake_minimum_required (VERSION 2.6.4)
include(FeatureSummary)

project (CPPQED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

##################################################
# Version
##################################################

# Versioning scheme

# MAJOR - the grand version of C++QED, as in v2
# MINOR - the milestone release, as in v2m9
# PATCH - should be increased on every snapshot package release, reset on milestone change

set(CPPQED_MAJOR_VERSION 2)
set(CPPQED_MINOR_VERSION 9)
set(CPPQED_PATCH_VERSION 0)
set(CPPQED_VERSION
  ${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION}.${CPPQED_PATCH_VERSION})


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

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
set(INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/CPPQED CACHE PATH
  "Installation directory for CMake files")
 
# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

##################################################
# Flags for release and debug mode
##################################################

# linker on MacOs doesn't like --no-undefined
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined -Wl,--as-needed")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-c++11-extensions -Wno-sign-compare -Wno-local-type-template-args")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG -std=c++11 -Wall -Wextra -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-deprecated -Wnon-virtual-dtor -Wno-ignored-qualifiers -Wno-c++11-extensions -Wno-sign-compare -Wno-local-type-template-args")


##################################################
# Compiler detection
##################################################

if (CMAKE_COMPILER_IS_GNUCXX)
  execute_process (
    COMMAND ${CMAKE_CXX_COMPILER} --version
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    RESULT_VARIABLE CMAKE_CXX_COMPILER_RESULT
    OUTPUT_VARIABLE CMAKE_CXX_COMPILER_OUTPUT
    ERROR_VARIABLE CMAKE_CXX_COMPILER_ERROR
  )
  if (CMAKE_CXX_COMPILER_RESULT EQUAL 0)
    string (REPLACE "\n" " " CMAKE_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_OUTPUT})
    set (CMAKE_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION} CACHE INTERNAL "C++ compiler version")
  endif (CMAKE_CXX_COMPILER_RESULT EQUAL 0)
endif (CMAKE_COMPILER_IS_GNUCXX)


##################################################
# Library detection
##################################################

find_package(GSL REQUIRED)
include_directories(include ${GSL_INCLUDE_DIRS})
set(THIRDPARTY_INCLUDE_DIRS ${THIRDPARTY_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

# The option SERIALIZATION will be disabled if any of the compoments required for serialization is missing
option(SERIALIZATION "Boost serialization" ON)
if(SERIALIZATION)
  set(ENABLE_SERIALIZATION 1)
else(SERIALIZATION)
  set(ENABLE_SERIALIZATION 0)
endif(SERIALIZATION)

# The option FLENS will be disabled if flens is missing
option(FLENS "FLENS support" ON)

# blitz
find_package(blitz REQUIRED)
if( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built with serialization support.")
else( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built without serialization support. Please configure blitz with --enable-serialization to enable optional binary statevector output.")
  set(ENABLE_SERIALIZATION 0)
endif( blitz_SERIALIZATION_FOUND )

include_directories( ${blitz_INCLUDE_DIR} )
set(THIRDPARTY_INCLUDE_DIRS ${THIRDPARTY_INCLUDE_DIRS} ${blitz_INCLUDE_DIR})


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
  set(ENABLE_SERIALIZATION 0)
endif( Boost_SERIALIZATION_FOUND )
include_directories( ${Boost_INCLUDE_DIRS} )
set(THIRDPARTY_INCLUDE_DIRS ${THIRDPARTY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

if(Boost_FOUND)
  if(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    add_definitions(-DBOOST_RESULT_OF_USE_TR1) 
    message(STATUS "added -DBOOST_RESULT_OF_USE_TR1" )
  endif()
endif()

# FLENS
find_package(flens)
if( flens_FOUND AND FLENS )
  include_directories(${flens_INCLUDE_DIRS})
  set(THIRDPARTY_INCLUDE_DIRS ${THIRDPARTY_INCLUDE_DIRS} ${flens_INCLUDE_DIRS})
  set(OPTLIBS ${OPTLIBS} ${flens_LIBRARIES})
  add_definitions(${flens_DEFINITIONS})
  set(ENABLE_FLENS 1)
else( flens_FOUND AND FLENS )
  message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  add_definitions( -DDO_NOT_USE_FLENS )
  set(ENABLE_FLENS 0)
endif( flens_FOUND AND FLENS )

# Check if serialization can be enabled, inform user
if( ENABLE_SERIALIZATION )
  message(STATUS "Support for binary statevector output enabled.")
  set(OPTLIBS ${OPTLIBS} ${Boost_SERIALIZATION_LIBRARY})
else( ENABLE_SERIALIZATION )
  add_definitions( -DDO_NOT_USE_BOOST_SERIALIZATION )
  message(STATUS "Optional support for binary statevector output disabled.")
endif( ENABLE_SERIALIZATION )

##################################################
# Compilation
##################################################

# build libutils.so
# add_subdirectory(utils)
# add utils/include to include path for all the other components
file(GLOB UTILS_SOURCE_FILES utils/src/*.cc)
include_directories(utils/include)

add_subdirectory(elements)
include_directories(${ELEMENTS_INCLUDE_DIRS})

# build libC++QEDcore.so
set(CORE_SOURCE_DIRS quantumdata quantumoperator structure quantumtrajectory)
foreach(d ${CORE_SOURCE_DIRS})
  aux_source_directory(${d} CORE_SOURCE_FILES)
  include_directories(${d})
endforeach(d)
set(CPPQEDLIB C++QED-${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION})
add_library(${CPPQEDLIB} SHARED ${CORE_SOURCE_FILES} ${ELEMENTS_SOURCE_FILES} ${UTILS_SOURCE_FILES})
target_link_libraries(${CPPQEDLIB} ${blitz_LIBRARY} ${OPTLIBS} ${GSL_LIBRARIES})
set_target_properties(${CPPQEDLIB}
		      PROPERTIES 
		      INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib
                      VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
                      SOVERSION ${CPPQED_ABI_MAJOR})


# build all the scripts
add_subdirectory(scripts)

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix)

##################################################
# Installation
##################################################

# The other libraries and scripts are installed in their own CMakeLists.txt
install(TARGETS ${CPPQEDLIB}
        EXPORT CPPQEDLibraryDepends
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
	COMPONENT shlib
)

# Install include files in their directory structure
set(ALL_INCLUDE_DIRS ${CORE_SOURCE_DIRS} ${ELEMENTS_INCLUDE_DIRS} utils/include)
foreach(d ${ALL_INCLUDE_DIRS})
  file(GLOB INC ${d}/*.h)
  file(GLOB DETAILS_INC ${d}/details/*.h)
  file(GLOB IMPL_INC ${d}/impl/*tcc)
  install(FILES ${INC} DESTINATION ${INSTALL_INCLUDE_DIR}/CPPQED/${d})
  install(FILES ${DETAILS_INC} DESTINATION ${INSTALL_INCLUDE_DIR}/CPPQED/${d}/details)
  install(FILES ${IMPL_INC} DESTINATION ${INSTALL_INCLUDE_DIR}/CPPQED/${d}/impl)
  file(RELATIVE_PATH NEW_INCLUDE_DIR ${INSTALL_CMAKE_DIR} ${INSTALL_INCLUDE_DIR}/CPPQED/${d})
  set(CONF_REL_INCLUDE_DIRS ${CONF_REL_INCLUDE_DIRS} ${NEW_INCLUDE_DIR})
endforeach(d)


# Add all targets to the build-tree export set
export(TARGETS ${CPPQEDLIB}
  FILE "${PROJECT_BINARY_DIR}/CPPQEDLibraryDepends.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
#export(PACKAGE CPPQED)
 
# Create a CPPQEDBuildTreeSettings.cmake file for the use from the build tree
configure_file(CPPQEDBuildTreeSettings.cmake.in
  "${PROJECT_BINARY_DIR}/CPPQEDBuildTreeSettings.cmake" @ONLY)
 
configure_file(CPPQEDConfig.cmake.in
  "${PROJECT_BINARY_DIR}/CPPQEDConfig.cmake" @ONLY)
configure_file(CPPQEDConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake" @ONLY)
 
# Install the CPPQEDConfig.cmake and CPPQEDConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/CPPQEDConfig.cmake"
  "${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake"
  "cmake/Modules/FindCBLAS.cmake"
  "cmake/Modules/Findflens.cmake"
  "cmake/Modules/Findblitz.cmake"
  "cmake/Modules/LibFindMacros.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT CPPQEDLibraryDepends DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

  
SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Raimar Sandner") #required
SET(CPACK_PACKAGE_CONTACT "vukics.andras@wigner.mta.hu")
SET(CPACK_PACKAGE_VERSION_MAJOR "${CPPQED_MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${CPPQED_MINOR_VERSION}")
#SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPPQED_VERSION}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPPQED_VERSION}")
 
INCLUDE(CPack)
  
##################################################
# Summary of enabled/disabled features
# (requires cmake >= 2.8.9)
##################################################

if(${CMAKE_VERSION} VERSION_GREATER "2.8.8")
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
  add_feature_info("FLENS" ENABLE_FLENS "compile framework with FLENS support.")
  add_feature_info(Serialization ENABLE_SERIALIZATION "needed for binary statevector output." )
  feature_summary( WHAT ALL )
endif(${CMAKE_VERSION} VERSION_GREATER "2.8.8")
