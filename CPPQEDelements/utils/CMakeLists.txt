# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDelements_utils OBJECT ModeFunction.cc)
add_library(CPPQED::CPPQEDelements_utils ALIAS CPPQEDelements_utils)

target_link_libraries(CPPQEDelements_utils PUBLIC CPPQED::CPPQEDcore)
target_include_directories(CPPQEDelements_utils PUBLIC .)
