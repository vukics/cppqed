# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_structure OBJECT DynamicsBase.cc ElementAveraged.cc)
add_library(CPPQED::CPPQEDcore_structure ALIAS CPPQEDcore_structure)

target_link_libraries(CPPQEDcore_structure PUBLIC CPPQED::CPPQEDcore_quantumdata)
target_include_directories(CPPQEDcore_structure PUBLIC .)
