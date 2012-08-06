cmake_minimum_required (VERSION 2.6)
include(FeatureSummary)

project (C++QED)

set(C++QED_CMAKE_MODULE_DIR ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Add packaged CMake scripts for backwards compatibility 
set(CMAKE_MODULE_PATH_INIT ${C++QED_CMAKE_MODULE_DIR})
file(GLOB BACKCOMP_VERSIONS RELATIVE ${C++QED_CMAKE_MODULE_DIR} ${C++QED_CMAKE_MODULE_DIR}/2.8.*)
foreach(BACKCOMP_VERSION ${BACKCOMP_VERSIONS})
  if (${CMAKE_VERSION} VERSION_LESS ${BACKCOMP_VERSION})
    set(CMAKE_MODULE_PATH_INIT ${CMAKE_MODULE_PATH_INIT} ${C++QED_CMAKE_MODULE_DIR}/${BACKCOMP_VERSION})
  endif()
endforeach()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_INIT} CACHE FILEPATH
  "Search path for CMake scripts and modules")

##################################################
# Flags for release and debug mode
##################################################

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG")

##################################################
# Library detection
##################################################

# This will be disabled if any of the compoments required for serialization is missing
option(WITH_SERIALIZATION "Boost serialization" ON)
if(WITH_SERIALIZATION)
  set(ENABLE_SERIALIATZION 1)
else(WITH_SERIALIZATION)
  set(ENABLE_SERIALIATZION 0)
endif(WITH_SERIALIZATION)

# blitz++

find_package(blitz++ REQUIRED)
if( blitz++_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built with serialization support.")
else( blitz++_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built without serialization support. Please configure blitz++ with --enable-serialization to enable optional binary statevector output.")
  set(ENABLE_SERIALIATZION 0)
endif( blitz++_SERIALIZATION_FOUND )

include_directories( ${blitz++_INCLUDE_DIR} )

set_package_properties(blitz++ PROPERTIES  
                                  DESCRIPTION "Multi-array library for C++."
                                  TYPE REQUIRED)
set_package_properties(blitz++_SERIALIZATION PROPERTIES  
                                  DESCRIPTION "Boost serialization support for blitz++ arrays."
                                  TYPE REQUIRED)

# Boost

find_package(Boost REQUIRED)
set_package_properties(Boost PROPERTIES  
                                  DESCRIPTION "Boost library."
                                  TYPE REQUIRED)
find_package(BoostSerialization)
if( BoostSerialization_FOUND )
  message(STATUS "Boost serialization library found.")
else( BoostSerialization_FOUND )
  message(STATUS "Boost serialization library not found.")
  set(ENABLE_SERIALIATZION 0)
endif( BoostSerialization_FOUND )
include_directories( ${Boost_INCLUDE_DIRS} )

# FLENS

find_package(flens)
if( flens_FOUND )
  include_directories(${flens_INCLUDE_DIR})
  set(OPTLIBS ${OPTLIBS} ${flens_LIBRARY})
else( flens_FOUND )
  message(STATUS "Flens library not found, optional flens support disabled.")
  add_definitions( -DDO_NOT_USE_FLENS )
endif( flens_FOUND )

# Check if serialization can be enabled, inform user

add_feature_info(Serialization ENABLE_SERIALIATZION "Support for binary statevector output." )

if( ENABLE_SERIALIZATION )
  #message(STATUS "Support for binary statevector output enabled.")
  add_definitions( -DUSE_BOOST_SERIALIZATION )
  set(OPTLIBS ${OPTLIBS} ${BoostSerialization_LIBRARY})
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

feature_summary( WHAT ALL )
