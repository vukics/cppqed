# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_high OBJECT QM_Picture.cc)
add_library(CPPQED::CPPQEDcore_high ALIAS CPPQEDcore_high)

target_link_libraries(CPPQEDcore_high PUBLIC CPPQED::CPPQEDcore_composites CPPQED::CPPQEDcore_quantumoperator CPPQED::CPPQEDcore_quantumtrajectory)
target_include_directories(CPPQEDcore_high PUBLIC .)
