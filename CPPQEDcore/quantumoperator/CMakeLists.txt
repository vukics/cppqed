# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_quantumoperator OBJECT Sigma.cc MultiDiagonal.cc)
add_library(CPPQED::CPPQEDcore_quantumoperator ALIAS CPPQEDcore_quantumoperator)

target_link_libraries(CPPQEDcore_quantumoperator PUBLIC CPPQED::CPPQEDcore_structure)
target_include_directories(CPPQEDcore_quantumoperator PUBLIC .)
