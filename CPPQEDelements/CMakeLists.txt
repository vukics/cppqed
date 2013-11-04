cmake_minimum_required (VERSION 2.8.9)

project(elements)

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

set(ELEMENTS_SOURCE_DIRS utils frees interactions)
set(frees_NEEDS utils)
set(interactions_NEEDS utils frees)

# This tells the elements_project macro that the project name 'elements' is ok
set(ORIGINAL_ELEMENTS_PROJECT 1)
elements_project()
