# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

add_library(CPPQEDutils SHARED
            BooleanNegatedProxy.cc CMatrix.cc Combinatorics.cc ComplexArrayExtensions.cc DrivenDampedHarmonicOscillator.cc Evolved.cc EvolvedGSL.cc FFT.cc FormDouble.cc
            FuzzyDouble.cc KeyPrinter.cc MathExtensions.cc Pars.cc ParsTrajectory.cc Randomized.cc Trajectory.cc Version.cc)

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC Blitz::Blitz)

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

# The serialization configuration of Blitz must be compatible with Boost_SERIALIZATION_FOUND
#if( Boost_SERIALIZATION_FOUND )

  #if(BZ_HAVE_BOOST_SERIALIZATION)
    #message(STATUS "Blitz++ built with serialization support.")
  #else()
    #message(FATAL_ERROR "Blitz++ built without serialization support. Please reconfigure blitz making sure it picks up Boost.Serialization.")
  #endif()

#endif( Boost_SERIALIZATION_FOUND )

#FLENS
##if( FLENS AND NOT FLENS_found )
  ##find_package(flens)
##endif( )

#if( NOT flens_FOUND )
  #message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  #set(DO_NOT_USE_FLENS ON)
#endif( )

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

#configure_file("${CPPQED_CMAKE_MODULE_PATH}/config.h.in" config.h)

#set(CPPQEDutils_LIB "CPPQEDutils-${CPPQED_ID}")
#file(GLOB CPPQEDutils_SOURCES . *.cc )

#add_library(${CPPQEDutils_LIB} SHARED ${CPPQEDutils_SOURCES})
#target_link_libraries(${CPPQEDutils_LIB} PUBLIC ${PUBLIC_LIBS} PRIVATE ${GSL_LIBRARIES})
#target_include_directories(${CPPQEDutils_LIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
#set_target_properties(${CPPQEDutils_LIB} PROPERTIES
#      PUBLIC_HEADER "${core_PUBLIC_HEADERS}"
#      INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
#      VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
#      SOVERSION ${CPPQED_ABI_MAJOR}
#)

#install(TARGETS ${CPPQEDutils_LIB} DESTINATION lib)

#file(GLOB CPPQEDutils_HEADERS . *.h )
#install(FILES ${CPPQEDutils_HEADERS} "${PROJECT_BINARY_DIR}/config.h" DESTINATION "include/CPPQEDutils-${CPPQED_ID}")
