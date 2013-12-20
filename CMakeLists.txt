cmake_minimum_required (VERSION 2.8.9)
project (CPPQED)

enable_testing()

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})
add_subdirectory(cpypyqed)

# Disabled CustomExample compilation by default for now, add it as a test later
# add_subdirectory(CustomElementsExample)
# set(CPPQEDelements_custom_example_DIR ${elements_custom_example_BINARY_DIR})
# add_subdirectory(CustomScriptsExample)

option(COMPILE_SCRIPTS "Compile the example scripts" On)
if(COMPILE_SCRIPTS)
  add_subdirectory(CPPQEDscripts)
endif()

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${CMAKE_INSTALL_DOCDIR}-${CPPQED_ID}
        PATTERN .git* EXCLUDE
        PATTERN .bzr EXCLUDE
        PATTERN *.kdev* EXCLUDE
        PATTERN build* EXCLUDE
)
add_subdirectory(Testing)