
#! \ingroup Main
#! \file
#! \brief Top level %CMake file for the C++QED elements component.
#!
#! The file structure is very simple, the main steps after finding the CPPQED core component are:
#! setting include dependency relations between sub-directories, calling elements_project() to
#! generate the library, and cppqed_documentation() to generate the documentation.


cmake_minimum_required (VERSION 2.8.9)

project(elements)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.100 REQUIRED)
include(${CPPQED_USE})
include_directories(${CPPQED_INCLUDE_DIRS})
############################################################

set(ELEMENTS_SOURCE_DIRS utils frees interactions)
set(frees_NEEDS utils)
set(interactions_NEEDS utils frees)

# This tells the elements_project macro that the project name 'elements' is ok
set(ORIGINAL_ELEMENTS_PROJECT 1)
elements_project()

##################################################
# Documentation
##################################################

if(CPPQED_MONOLITHIC)
  cppqed_documentation("elements_" ${core_DOXYGEN_TAG} core_doc)
endif()
