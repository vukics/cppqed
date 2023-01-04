# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

# ### Library detection
find_package(Boost ${BOOST_MINIMAL} REQUIRED COMPONENTS iostreams json serialization)

find_package(Eigen3 3.3 REQUIRED)

add_subdirectory(thirdParty)

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            #BooleanNegatedProxy.cc CMatrix.cc ComplexArrayExtensions.cc DrivenDampedHarmonicOscillator.cc FFT.cc FormDouble.cc FuzzyDouble.cc KeyPrinter.cc MathExtensions.cc
            Pars.cc Trajectory.cc FormDouble.cc KeyPrinter.cc
            #Random.cc 
            Version.cc
            )

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC ThirdParty Boost::iostreams Boost::json Boost::serialization Eigen3::Eigen)

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

