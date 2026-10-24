# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDelements_subsystems OBJECT
            Mode.cc
#            BichromaticMode.cc HomodynedMode.cc ParsMode.cc
#            ParsParticle.cc Particle.cc DrivenTwoLevelAtom.cc
            Qbit.cc Spin.cc # TimeIndependentMatrixHamiltonian.cc
            )
add_library(CPPQED::CPPQEDelements_subsystems ALIAS CPPQEDelements_subsystems)

target_link_libraries(CPPQEDelements_subsystems PUBLIC CPPQED::CPPQEDcore)
target_include_directories(CPPQEDelements_subsystems PUBLIC .)
