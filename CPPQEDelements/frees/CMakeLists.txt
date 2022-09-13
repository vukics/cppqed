# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDelements_frees OBJECT
            BichromaticMode.cc HomodynedMode.cc Mode.cc ParsMode.cc
            ParsParticle.cc Particle.cc DrivenTwoLevelAtom.cc
            Qbit.cc Spin.cc TimeIndependentMatrixHamiltonian.cc)
add_library(CPPQED::CPPQEDelements_frees ALIAS CPPQEDelements_frees)

target_link_libraries(CPPQEDelements_frees PUBLIC CPPQED::CPPQEDelements_utils)
target_include_directories(CPPQEDelements_frees PUBLIC .)
