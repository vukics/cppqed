# Copyright András Vukics 2020–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDutils")

# ### Library detection
find_package(Boost ${BOOST_MINIMAL} REQUIRED COMPONENTS iostreams json serialization)

find_package(Eigen3 3.4 REQUIRED)

add_subdirectory(thirdParty)

add_library(CPPQEDutils SHARED
            $<TARGET_OBJECTS:cppqedGitSHA1>
            DrivenDampedHarmonicOscillator.cc Pars.cc Trajectory.cc Version.cc
            )

add_library(CPPQED::CPPQEDutils ALIAS CPPQEDutils)

target_link_libraries(CPPQEDutils PUBLIC ThirdParty Boost::json Boost::serialization Eigen3::Eigen)

install(TARGETS CPPQEDutils DESTINATION lib)