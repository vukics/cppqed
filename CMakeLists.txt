cmake_minimum_required (VERSION 2.8.9)

project(CPPQEDelements)

############################################################
# This section is essential to find the required libraries 
# and include files
if(NOT DEFINED CPPQED_MONOLITHIC)
  find_package(CPPQED 2.9 REQUIRED)
else(NOT DEFINED CPPQED_MONOLITHIC)
  include(${CPPQED_BINARY_DIR}/CPPQEDcore/CPPQEDConfig.cmake)
endif(NOT DEFINED CPPQED_MONOLITHIC)
include(${CPPQED_USE})
CPPQED_SETUP()
############################################################

add_subdirectory(src)