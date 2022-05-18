# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

# ### Library detection
find_package(GSL 2.5 REQUIRED)

# Boost
find_package(Boost ${BOOST_MINIMAL} REQUIRED COMPONENTS iostreams json)

find_package(Eigen3 3.3 REQUIRED)

add_subdirectory(thirdPartyRandom)

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            #BooleanNegatedProxy.cc CMatrix.cc ComplexArrayExtensions.cc DrivenDampedHarmonicOscillator.cc FFT.cc FormDouble.cc FuzzyDouble.cc KeyPrinter.cc MathExtensions.cc Pars.cc Random.cc Trajectory.cc Version.cc
            )

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC ThirdPartyRandom Boost::iostreams Boost::json GSL::gsl Eigen3::Eigen)

target_compile_definitions(CPPQEDutils PUBLIC CPPQED_HAS_GSL)

# Project options
# Switch for boost bzip2 support
option(COMPRESSION "Boost bzip2 support" ON)


# Check if compression can be enabled, inform user
if(COMPRESSION)
  message(STATUS "Support for binary statevector compression enabled.")
else()
  message(STATUS "Optional support for binary statevector compression disabled.")
  target_compile_definitions(CPPQEDutils PUBLIC DO_NOT_USE_BOOST_COMPRESSION)
endif()

add_subdirectory(examples)

if(ENABLE_TESTING)
  add_subdirectory(Testing)
endif()

add_subdirectory(sandbox)

install(TARGETS CPPQEDutils DESTINATION lib)

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
