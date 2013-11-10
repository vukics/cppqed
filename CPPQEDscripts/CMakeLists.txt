cmake_minimum_required (VERSION 2.8.9)

project(scripts)

############################################################
# This section is essential to find the required libraries 
# and include files
if(NOT DEFINED CPPQED_MONOLITHIC)
  find_package(CPPQED 2.99 REQUIRED)
else(NOT DEFINED CPPQED_MONOLITHIC)
  include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
endif(NOT DEFINED CPPQED_MONOLITHIC)
include(${CPPQED_USE})
############################################################

set(NEED_FLENS GeneralDickeImaginaryEvolution)
set(EXCLUDE_SCRIPTS GeneralDickeImaginaryEvolution Ca40InCavity Ca40InCavityWithPhonons CavityRaman MLJC NX_coupledModesElim PumpedLossyModeRegression Raman Spin)
scripts_project()

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix Ring)
