# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDcore")

add_subdirectory(quantumdata)

add_subdirectory(structure)

add_subdirectory(quantumoperator)

add_subdirectory(quantumtrajectory)

add_library(CPPQEDcore SHARED
            $<TARGET_OBJECTS:CPPQED::CPPQEDcore_quantumoperator>)
add_library(CPPQED::CPPQEDcore ALIAS CPPQEDcore)

target_link_libraries(CPPQEDcore PUBLIC CPPQED::CPPQEDcore_quantumoperator CPPQED::CPPQEDcore_quantumtrajectory)

install(TARGETS CPPQEDcore DESTINATION lib)