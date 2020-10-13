# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDscripts")

set(FEWER_SCRIPTS_TARGETS
    PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED
    QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix SingleParticle Ring Ring_Evolved Raman CavityRaman)

foreach(S ${FEWER_SCRIPTS_TARGETS})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDelements)
endforeach()
    
# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts ${FEWER_SCRIPTS_TARGETS})

#set(NEED_FLENS GeneralDickeImaginaryEvolution)
# If something needs to be excluded, that’s how to do it:
# set(EXCLUDE_SCRIPTS GeneralDickeImaginaryEvolution NX_coupledModesElim PumpedLossyModeRegression Spin)
# If the target should exist, but the script should not be compiled as part of ALL, add it to the following
# variable:
# set(EXCLUDE_FROM_ALL_SCRIPTS 6qbits)

#set(EXCLUDE_SCRIPTS Ca40InCavityWithPhonons)
#scripts_project()

