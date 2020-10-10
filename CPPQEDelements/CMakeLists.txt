# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDelements")

add_subdirectory(utils)

add_subdirectory(frees)

add_subdirectory(interactions)

add_library(CPPQEDelements INTERFACE)
add_library(CPPQED::CPPQEDelements ALIAS CPPQEDelements)

target_link_libraries(CPPQEDelements INTERFACE CPPQED::CPPQEDelements_interactions)
