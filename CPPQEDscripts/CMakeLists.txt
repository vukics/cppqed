cmake_minimum_required (VERSION 2.6.4)

project(CPPQEDscripts)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.8 REQUIRED)
include(CPPQEDCommonChecks)
CPPQED_CXX_FLAGS()
include_directories(${CPPQED_INCLUDE_DIRS})
include_directories(SYSTEM ${CPPQED_THIRDPARTY_INCLUDE_DIRS})
add_definitions(${CPPQED_DEFINITIONS})
############################################################

message(STATUS FLENS: ${CPPQED_HAS_FLENS})
message(STATUS Serialization: ${CPPQED_HAS_SERIALIZATION})

add_subdirectory(src)