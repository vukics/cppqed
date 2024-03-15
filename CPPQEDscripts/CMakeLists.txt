# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDscripts")

# now this is the only collective scripts target, compiles everything needed for the testsuite
# think over: fewer_scripts should contain such scripts that collectively exemplify all the important features of the framework,
# so that if they compile then at least the compile-time features are ok
set(FEWER_SCRIPTS_TARGETS
    DTLA_Evolved DTLA_C++QED DrivenDissipativeQbit DTLA_EvolvedHL DrivenDissipativeMode_Evolved DrivenDissipativeMode
    QbitMode QbitMode_Evolved QbitMode_Matrix SingleParticle 1particle1mode Ring Ring_Evolved Raman CavityRaman)
    
foreach(S ${FEWER_SCRIPTS_TARGETS})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDelements)
  install(TARGETS ${S} DESTINATION bin)
endforeach()

add_executable(Spin Spin.cc)
target_link_libraries(Spin CPPQED::CPPQEDelements)

# add target "fewer_scripts"
add_custom_target(fewer_scripts DEPENDS ${FEWER_SCRIPTS_TARGETS})

#set(NEED_FLENS GeneralDickeImaginaryEvolution)
# If something needs to be excluded, that’s how to do it:
# set(EXCLUDE_SCRIPTS GeneralDickeImaginaryEvolution NX_coupledModesElim DrivenDissipativeModeRegression Spin)
# If the target should exist, but the script should not be compiled as part of ALL, add it to the following
# variable:
# set(EXCLUDE_FROM_ALL_SCRIPTS 6qbits)

#set(EXCLUDE_SCRIPTS Ca40InCavityWithPhonons)
#scripts_project()

