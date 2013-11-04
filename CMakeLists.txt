cmake_minimum_required (VERSION 2.8.9)

project (CPPQED)

set(CPPQED_MONOLITHIC 1)

add_subdirectory(CPPQEDcore)
add_subdirectory(CPPQEDelements)
add_subdirectory(CPPQEDscripts)