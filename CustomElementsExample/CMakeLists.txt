cmake_minimum_required (VERSION 2.8.9)

# The project name must be unique and cannot be "elements".
# This is ne name by which we refer to the library in scripts projects.
project(elements_custom_example)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.100 REQUIRED)
include(${CPPQED_USE})
############################################################

# Define all source directories and their include dependencies, e.g.:
#
# set(ELEMENTS_SOURCE_DIRS utils frees interactions)
# set(frees_NEEDS utils)
# set(interactions_NEEDS utils frees)
#
# Put a CMakeLists.txt file in each subdirectory with the single line
# create_object_target()


# In this simple example, there is only one subdirectory and hence no dependencies
set(ELEMENTS_SOURCE_DIRS interactions)

# all the magic happens here
elements_project()
