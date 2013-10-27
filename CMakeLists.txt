cmake_minimum_required (VERSION 2.8.9)

project(CPPQEDelements)

include(GNUInstallDirs)

############################################################
# This section is essential to find the required libraries 
# and include files
if(NOT DEFINED CPPQED_MONOLITHIC)
  find_package(CPPQED 2.9 REQUIRED)
else(NOT DEFINED CPPQED_MONOLITHIC)
  include(${CPPQED_BINARY_DIR}/CPPQEDcore/CPPQEDConfig.cmake)
endif(NOT DEFINED CPPQED_MONOLITHIC)
include(${CPPQED_USE})
CPPQED_CXX_FLAGS()
include_directories(${CPPQED_INCLUDE_DIRS})
include_directories(SYSTEM ${CPPQED_THIRDPARTY_INCLUDE_DIRS})
############################################################


add_subdirectory(src)