cmake_minimum_required (VERSION 2.8.9)

project(elements_custom_example)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.9 REQUIRED)
include(${CPPQED_USE})
############################################################

set(ELEMENTS_SOURCE_DIRS interactions)

elements_project()
