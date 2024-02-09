# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_quantumtrajectory OBJECT EvolutionMethod.cc)
add_library(CPPQED::CPPQEDcore_quantumtrajectory ALIAS CPPQEDcore_quantumtrajectory)

target_link_libraries(CPPQEDcore_quantumtrajectory PUBLIC CPPQED::CPPQEDcore_structure)
target_include_directories(CPPQEDcore_quantumtrajectory PUBLIC .)
