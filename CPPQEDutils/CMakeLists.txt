# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            BooleanNegatedProxy.cc CMatrix.cc Combinatorics.cc ComplexArrayExtensions.cc Evolved.cc EvolvedGSL.cc FFT.cc FormDouble.cc
            FuzzyDouble.cc KeyPrinter.cc MathExtensions.cc Pars.cc ParsTrajectory.cc Randomized.cc Trajectory.cc Version.cc)
# DrivenDampedHarmonicOscillator.cc should have an FLENS-dependent logic

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

# ### Library detection
find_package(GSL 2.5 REQUIRED)

# Boost
find_package(Boost 1.71 REQUIRED COMPONENTS iostreams)

target_link_libraries(CPPQEDutils PUBLIC Blitz::Blitz Boost::iostreams PRIVATE GSL::gsl)

# Project options
# Switch for boost bzip2 support
option(COMPRESSION "Boost bzip2 support" ON)

#FLENS
##if( FLENS AND NOT FLENS_found )
  ##find_package(flens)
##endif( )

#if( NOT flens_FOUND )
  #message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  #set(DO_NOT_USE_FLENS ON)
#endif( )

# include_directories(SYSTEM ${flens_INCLUDE_DIRS})

# Check if compression can be enabled, inform user
if(COMPRESSION)
  message(STATUS "Support for binary statevector compression enabled.")
else()
  set(DO_NOT_USE_BOOST_COMPRESSION ON)
  message(STATUS "Optional support for binary statevector compression disabled.")
endif()

if(ENABLE_TESTING)
  add_subdirectory(Testing)
endif()

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
