# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

# ### Library detection
find_package(Boost ${BOOST_MINIMAL} REQUIRED COMPONENTS iostreams serialization)

find_package(Eigen3 3.3 REQUIRED)

add_subdirectory(thirdParty)

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            DrivenDampedHarmonicOscillator.cc Pars.cc Trajectory.cc Version.cc
            )

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC ThirdParty Boost::iostreams Boost::serialization Eigen3::Eigen)

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

install(TARGETS CPPQEDutils DESTINATION lib)

