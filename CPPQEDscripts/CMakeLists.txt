cmake_minimum_required (VERSION 2.8.9)

project(scripts)

############################################################
# This section is essential to find the required libraries 
# and include files
if(NOT DEFINED CPPQED_MONOLITHIC)
  find_package(CPPQED 2.99 REQUIRED)
else(NOT DEFINED CPPQED_MONOLITHIC)
  include(${CPPQED_BINARY_DIR}/CPPQEDcore/CPPQEDConfig.cmake)
endif(NOT DEFINED CPPQED_MONOLITHIC)
include(${CPPQED_USE})
############################################################

set(NEED_FLENS GeneralDickeImaginaryEvolution)
set(EXCLUDE_SCRIPTS Ca40InCavity Ca40InCavityWithPhonons CavityRaman MLJC NX_coupledModesElim PumpedLossyModeRegression Raman Spin)
scripts_project()