cmake_minimum_required (VERSION 2.8.9)
project (CPPQED)

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

add_subdirectory(CPPQEDcore)
add_subdirectory(CPPQEDelements)

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