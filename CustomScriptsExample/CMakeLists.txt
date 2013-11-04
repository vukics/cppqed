cmake_minimum_required (VERSION 2.6.4)

project(scripts_custom_example)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.9 REQUIRED)
include(${CPPQED_USE})
############################################################

set(scripts_custom_example_ELEMENTS elements_custom_example)
scripts_project()
