cmake_minimum_required (VERSION 2.8.9)

project(scripts_custom_example)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.10 REQUIRED)
include(${CPPQED_USE})
############################################################

# All *.cc files in the toplevel directory will be picked up.

# scripts_project(elements_project1 elements_project2 ...)
# will link to elements_project1, elements_project2 etc
# Note that CPPQEDelements is used automatically
scripts_project(elements_custom_example)
