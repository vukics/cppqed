cmake_minimum_required (VERSION 2.8.9)

project(scripts)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.99 REQUIRED)
include(${CPPQED_USE})
############################################################

set(NEED_FLENS GeneralDickeImaginaryEvolution)
set(EXCLUDE_SCRIPTS GeneralDickeImaginaryEvolution Ca40InCavity Ca40InCavityWithPhonons CavityRaman MLJC NX_coupledModesElim PumpedLossyModeRegression Raman Spin)
# Exclude 4qbits scripts form the "all" target for now, as they take much memory to compile, which can
# cause excessive swapping when building things in parallel. Also, they are suspected to cause random g++ crashes 
set(EXCLUDE_FROM_ALL_SCRIPTS 4qbits 4qbitsWF 3qbits 3qbitsPostprocess)
scripts_project()

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix Ring)
