
#! \ingroup Main
#! \file
#! \brief Top level %CMake file for the C++QED scripts component.
#!
#! The file structure is very simple, all the work is done by the scripts_project() macro.


cmake_minimum_required (VERSION 2.8.9)

project(scripts)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.100 REQUIRED)
include(${CPPQED_USE})
############################################################

set(NEED_FLENS GeneralDickeImaginaryEvolution)
# If something needs to be excluded, that’s how to do it:
# set(EXCLUDE_SCRIPTS GeneralDickeImaginaryEvolution NX_coupledModesElim PumpedLossyModeRegression Spin)
# If the target should exist, but the script should not be compiled as part of ALL, add it to the following
# variable:
# set(EXCLUDE_FROM_ALL_SCRIPTS 6qbits)

set(EXCLUDE_SCRIPTS Ca40InCavityWithPhonons CavityRaman)
scripts_project()

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix SingleParticle Ring Ring_Evolved)
