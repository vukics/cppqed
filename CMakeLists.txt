cmake_minimum_required (VERSION 2.6.4)
include(FeatureSummary)

project (C++QED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

##################################################
# Flags for release and debug mode
##################################################

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG -Wall -Wextra -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnon-virtual-dtor")

##################################################
# Library detection
##################################################

# This will be disabled if any of the compoments required for serialization is missing
option(SERIALIZATION "Boost serialization" ON)
if(SERIALIZATION)
  set(ENABLE_SERIALIZATION 1)
else(SERIALIZATION)
  set(ENABLE_SERIALIZATION 0)
endif(SERIALIZATION)

option(FLENS "FLENS support" ON)

# blitz++

find_package(blitz++ REQUIRED)
if( blitz++_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built with serialization support.")
else( blitz++_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built without serialization support. Please configure blitz++ with --enable-serialization to enable optional binary statevector output.")
  set(ENABLE_SERIALIZATION 0)
endif( blitz++_SERIALIZATION_FOUND )

include_directories( ${blitz++_INCLUDE_DIR} )


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

# FLENS

find_package(flens)
if( flens_FOUND AND FLENS )
  include_directories(${flens_INCLUDE_DIR})
  set(OPTLIBS ${OPTLIBS} ${flens_LIBRARIES})
  add_definitions(${flens_DEFINITIONS})
else( flens_FOUND AND FLENS )
  message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  add_definitions( -DDO_NOT_USE_FLENS )
endif( flens_FOUND AND FLENS )

# Check if serialization can be enabled, inform user

if( ENABLE_SERIALIZATION )
  message(STATUS "Support for binary statevector output enabled.")
  add_definitions( -DUSE_BOOST_SERIALIZATION )
  set(OPTLIBS ${OPTLIBS} ${Boost_SERIALIZATION_LIBRARY})
else( ENABLE_SERIALIZATION )
  message(STATUS "Optional support for binary statevector output disabled.")
endif( ENABLE_SERIALIZATION )

##################################################
# Compilation
##################################################

# build libutils.so
add_subdirectory(utils)
# add utils/include to include path for all the other components
include_directories(${C++QED_SOURCE_DIR}/utils/include)

# build libC++QED.so
set(CORE_SOURCE_DIRS quantumdata quantumoperator structure quantumtrajectory)
foreach(d ${CORE_SOURCE_DIRS})
  aux_source_directory(${d} CORE_SOURCE_FILES)
  include_directories(${d})
endforeach(d)
add_library(C++QED SHARED ${CORE_SOURCE_FILES})
target_link_libraries( C++QED utils ${blitz++_LIBRARY} ${OPTLIBS} )

# build libel.so
add_subdirectory(elements)

# build all the scripts
add_subdirectory(scripts)

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix)


##################################################
# Summary of enabled/disabled features
# (requires cmake >= )
##################################################

if(${CMAKE_VERSION} VERSION_GREATER "2.8.7")
  set_package_properties(PkgConfig PROPERTIES URL "http://pkgconfig.freedesktop.org/wiki"
                                 DESCRIPTION "Package config system that manages compile/link flags"
                                 TYPE OPTIONAL
                                 PURPOSE "Assists cmake in finding libraries.")
  set_package_properties(blitz++ PROPERTIES  URL "http://sf.net/projects/blitz/"
                                 DESCRIPTION "High-performance C++ vector mathematics library"
                                 TYPE REQUIRED
                                 PURPOSE "Multi-Array implementation used in the framework.")
  set_package_properties(flens PROPERTIES URL "http://www.mathematik.uni-ulm.de/~lehn/FLENS/"
                                 DESCRIPTION "Flexible Library for Efficient Numerical Solutions."
                                 TYPE OPTIONAL
                                 PURPOSE "With FLENS the partial transpose of the density operator can be calculated.")
  set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                 DESCRIPTION "Boost C++ Library."
                                 TYPE REQUIRED
                                 PURPOSE "Advanced template metaprogramming and preprocessor algorithms.")
  set_package_properties(GSL PROPERTIES URL "http://www.gnu.org/software/gsl/"
                                 DESCRIPTION "GNU Scientific Library"
                                 TYPE REQUIRED
                                 PURPOSE "Used as implementation of an ODE solver and random number generator.")
  set_package_properties(CBLAS PROPERTIES URL "http://www.gnu.org/software/gsl/"
                                 DESCRIPTION "C library for Basic Linear Algebra Subprograms"
                                 TYPE REQUIRED)
  add_feature_info(Serialization ENABLE_SERIALIZATION "needed for binary statevector output." )
  feature_summary( WHAT ALL )
endif(${CMAKE_VERSION} VERSION_GREATER "2.8.7")
