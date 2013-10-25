cmake_minimum_required (VERSION 2.6.4)

project(CPPQEDelements)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.9 REQUIRED)
include(${CPPQED_USE})
CPPQED_CXX_FLAGS()
message(STATUS ${CPPQED_INCLUDE_DIRS})
include_directories(${CPPQED_INCLUDE_DIRS})
include_directories(SYSTEM ${CPPQED_THIRDPARTY_INCLUDE_DIRS})
############################################################
add_subdirectory(src)