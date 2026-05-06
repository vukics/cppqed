# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_quantumtrajectory INTERFACE)
add_library(CPPQED::CPPQEDcore_quantumtrajectory ALIAS CPPQEDcore_quantumtrajectory)

target_link_libraries(CPPQEDcore_quantumtrajectory INTERFACE CPPQED::CPPQEDcore_structure)
target_include_directories(CPPQEDcore_quantumtrajectory INTERFACE .)