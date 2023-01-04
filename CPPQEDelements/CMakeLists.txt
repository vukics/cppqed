# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDelements")

add_subdirectory(utils)

add_subdirectory(frees)

add_subdirectory(interactions)

add_library(CPPQEDelements SHARED 
            $<TARGET_OBJECTS:CPPQED::CPPQEDelements_utils> $<TARGET_OBJECTS:CPPQED::CPPQEDelements_frees> $<TARGET_OBJECTS:CPPQED::CPPQEDelements_interactions>)
add_library(CPPQED::CPPQEDelements ALIAS CPPQEDelements)

target_link_libraries(CPPQEDelements PUBLIC CPPQED::CPPQEDelements_interactions)

install(TARGETS CPPQEDelements DESTINATION lib)
