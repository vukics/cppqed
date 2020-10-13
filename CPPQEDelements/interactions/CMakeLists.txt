# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDelements_interactions OBJECT
            CoupledModes.cc GeneralDicke.cc JaynesCummings.cc ModeCorrelations.cc
            ParsParticleCavity.cc ParticleCavity.cc ParticleTwoModes.cc QbitModeCorrelations.cc)
add_library(CPPQED::CPPQEDelements_interactions ALIAS CPPQEDelements_interactions)

target_link_libraries(CPPQEDelements_interactions PUBLIC CPPQED::CPPQEDelements_frees)
target_include_directories(CPPQEDelements_interactions PUBLIC .)
