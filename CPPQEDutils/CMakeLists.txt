# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

# If BUNDLED_BLITZ is ON, then we know in advance where blitz *will be* situated, so we can skip FindBlitz, and simply set the blitz variables by hand
# thereupon, we have a unified treatment available, independently of the value of BUNDLED_BLITZ
# so, we can confine BUNDLED_BLITZ to the highest level only (main CMakeLists.txt)

# Ditto for FLENS

project(cppqedutils)

message(STATUS "Configuring CPPQEDutils")

set(CPPQED_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
list(APPEND CMAKE_MODULE_PATH ${CPPQED_CMAKE_MODULE_PATH})

#! \file
#! <!--#########################################################-->
#! ### Library detection
#! <!--#########################################################-->
#!
#! In this section we look for required and optional dependencies.

find_package(GSL 2.5 REQUIRED)

# Boost
find_package(Boost 1.71 REQUIRED COMPONENTS iostreams OPTIONAL_COMPONENTS serialization)

#! \name Project options
#! @{

#! Switch for boost bzip2 support
option(COMPRESSION "Boost bzip2 support" ON)

#! @}

# Blitz
if( NOT blitz_FOUND )
  find_package(blitz REQUIRED)
endif( NOT blitz_FOUND )

list(APPEND PUBLIC_LIBS ${Boost_IOSTREAMS_LIBRARY} ${blitz_LIBRARIES})

# The serialization configuration of Blitz must be compatible with Boost_SERIALIZATION_FOUND
if( Boost_SERIALIZATION_FOUND )

  if( blitz_SERIALIZATION_FOUND )
    message(STATUS "Blitz++ built with serialization support.")
  else( blitz_SERIALIZATION_FOUND )
    message(FATAL_ERROR "Blitz++ built without serialization support. Please reconfigure blitz making sure it picks up Boost.Serialization.")
  endif( blitz_SERIALIZATION_FOUND )

endif( Boost_SERIALIZATION_FOUND )

#FLENS
if( FLENS AND NOT FLENS_found )
  find_package(flens)
endif( )

if( NOT flens_FOUND )
  message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  set(DO_NOT_USE_FLENS ON)
endif( )

include_directories(SYSTEM ${blitz_INCLUDE_DIRS} ${flens_INCLUDE_DIRS})

# Check if serialization can be enabled, inform user
if( Boost_SERIALIZATION_FOUND )
  message(STATUS "Support for binary statevector output enabled.")
  list(APPEND PUBLIC_LIBS ${Boost_SERIALIZATION_LIBRARY})
else( Boost_SERIALIZATION_FOUND )
  set(DO_NOT_USE_BOOST_SERIALIZATION ON)
  message(STATUS "Optional support for binary statevector output disabled.")
endif( Boost_SERIALIZATION_FOUND )

# Check if compression can be enabled, inform user
if( COMPRESSION )
  message(STATUS "Support for binary statevector compression enabled.")
else( COMPRESSION )
  set(DO_NOT_USE_BOOST_COMPRESSION ON)
  message(STATUS "Optional support for binary statevector compression disabled.")
endif( COMPRESSION )

configure_file("${CPPQED_CMAKE_MODULE_PATH}/config.h.in" config.h)

file(GLOB CPPQEDutils_sources . *.cc )
list(REMOVE_ITEM CPPQEDutils_sources "Version.cc")
add_library(CPPQEDutils SHARED ${CPPQEDutils_sources})
target_link_libraries(CPPQEDutils PUBLIC ${PUBLIC_LIBS} PRIVATE ${GSL_LIBRARIES})
target_include_directories(CPPQEDutils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
# create_object_target()
