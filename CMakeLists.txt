cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project (C++QED)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG")

find_package(blitz++ REQUIRED)
find_package(Boost REQUIRED)
find_package(Boost QUIET COMPONENTS serialization)
find_package(flens)

set(ENABLE_SERIALIATZION 1)
if( BLITZ_SERIALIZATION )
  message("Blitz++ built with serialization support.")
else( BLITZ_SERIALIZATION )
  message("Blitz++ built without serialization support. Please configure blitz++ with --enable-serialization to enable optional binary statevector output.")
  set(ENABLE_SERIALIATZION 0)
endif( BLITZ_SERIALIZATION )

if( Boost_SERIALIZATION_FOUND )
  message("Boost serialization library found.")
else( Boost_SERIALIZATION_FOUND )
  message("Boost serialization library not found.")
  set(ENABLE_SERIALIATZION 0)
endif( Boost_SERIALIZATION_FOUND )

if( ENABLE_SERIALIATZION )
  message("Support for binary statevector output enabled.")
  add_definitions( -DUSE_BOOST_SERIALIZATION )
else( ENABLE_SERIALIATZION )
  message("Optional support for binary statevector output disabled.")
endif( ENABLE_SERIALIATZION )

if( flens_FOUND )
  include_directories(${flens_INCLUDE_DIR})
else( flens_FOUND )
  message("Flens library not found, optional flens support disabled.")
  add_definitions( -DDO_NOT_USE_FLENS )
endif( flens_FOUND )

set(INCLUDEDIRS ${blitz++_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${flens_INCLUDE_DIR})

add_subdirectory(utils)
include_directories(${C++QED_SOURCE_DIR}/utils/include)

set(CORE_SOURCE_DIRS quantumdata quantumoperator structure quantumtrajectory)

foreach(d ${CORE_SOURCE_DIRS})
  aux_source_directory(${d} CORE_SOURCE_FILES)
  include_directories(${d})
endforeach(d)

add_library(C++QED SHARED ${CORE_SOURCE_FILES})
target_link_libraries(C++QED utils ${blitz++_LIBRARY} ${flens_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} )

add_subdirectory(elements)
add_subdirectory(scripts)

add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix)
