# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_composites OBJECT BinarySystem.cc)
add_library(CPPQED::CPPQEDcore_composites ALIAS CPPQEDcore_composites)

target_link_libraries(CPPQEDcore_composites PUBLIC CPPQED::CPPQEDcore_structure)
target_include_directories(CPPQEDcore_composites PUBLIC .)
