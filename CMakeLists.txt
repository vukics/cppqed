cmake_minimum_required (VERSION 2.8.9)
project (CPPQED)

enable_testing()

include(GNUInstallDirs)

set(CPPQED_MONOLITHIC 1)

set(DOC_INSTALL_DIR ${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID})

add_subdirectory(CPPQEDcore)
set(CPPQED_DIR ${core_BINARY_DIR})
add_subdirectory(CPPQEDelements)
set(CPPQEDelements_DIR ${elements_BINARY_DIR})

find_package(Boost QUIET COMPONENTS python)
if(Boost_PYTHON_FOUND)
  add_subdirectory(cpypyqed)
endif()

# Disabled CustomExample compilation by default for now, add it as a test later
# add_subdirectory(CustomElementsExample)
# set(CPPQEDelements_custom_example_DIR ${elements_custom_example_BINARY_DIR})
# add_subdirectory(CustomScriptsExample)

option(COMPILE_SCRIPTS "Compile the example scripts" On)
if(COMPILE_SCRIPTS)
  add_subdirectory(CPPQEDscripts)
endif()

install(DIRECTORY CustomElementsExample CustomScriptsExample
        DESTINATION ${DOC_INSTALL_DIR}
        PATTERN .git* EXCLUDE
        PATTERN .bzr EXCLUDE
        PATTERN *.kdev* EXCLUDE
        PATTERN build* EXCLUDE
)
add_subdirectory(Testing)

##################################################
# Documentation
##################################################

find_package(Doxygen)

if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  set(CONF_DOC_DIR ${CMAKE_BINARY_DIR}/doc)
  file(MAKE_DIRECTORY ${CONF_DOC_DIR})
  file(RELATIVE_PATH CPPQED_RELATIVE_DOXYGEN_TAG ${CONF_DOC_DIR} ${CONF_DOC_DIR}/core/core.tag)
  get_filename_component(CPPQED_RELATIVE_DOXYGEN_DIR ${CPPQED_RELATIVE_DOXYGEN_TAG} DIRECTORY)
  set(CPPQED_RELATIVE_DOXYGEN_DIR ../${CPPQED_RELATIVE_DOXYGEN_DIR})
  configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${CONF_DOC_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CONF_DOC_DIR}/Doxyfile
    WORKING_DIRECTORY ${CONF_DOC_DIR}
    COMMENT "Generating APIs documentation with Doxygen"
    DEPENDS core_doc elements_doc
  )
  install(DIRECTORY ${CONF_DOC_DIR}/html
        DESTINATION ${DOC_INSTALL_DIR}
        OPTIONAL
  )
endif()