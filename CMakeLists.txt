cmake_minimum_required (VERSION 2.8.9)

project(elements)

############################################################
# This section is essential to find the required libraries 
# and include files
find_package(CPPQED 2.99 REQUIRED)
include(${CPPQED_USE})
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

find_package(Doxygen)

if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  if(CPPQED_MONOLITHIC)
    set(DOXYGEN_HEADER_FILE ${cppqed_SOURCE_DIR}/doc/header.html)
    set(DOXYGEN_CSS_FILE ${cppqed_SOURCE_DIR}/doc/stylesheet.css)
  endif()
  set(CONF_DOC_DIR ${CMAKE_BINARY_DIR}/doc/${PROJECT_NAME})
  file(MAKE_DIRECTORY ${CONF_DOC_DIR})
  file(RELATIVE_PATH CPPQED_RELATIVE_DOXYGEN_TAG ${CONF_DOC_DIR} ${CPPQED_DOXYGEN_TAG})
  get_filename_component(CPPQED_RELATIVE_DOXYGEN_DIR ${CPPQED_RELATIVE_DOXYGEN_TAG} PATH)
  configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in ${CONF_DOC_DIR}/Doxyfile @ONLY)
  add_custom_target(${PROJECT_NAME}_doc ${DOXYGEN_EXECUTABLE} ${CONF_DOC_DIR}/Doxyfile
    WORKING_DIRECTORY ${CONF_DOC_DIR}
    COMMENT "Generating API documentation with Doxygen"
  )
  if(CPPQED_MONOLITHIC)
    add_dependencies(${PROJECT_NAME}_doc core_doc)
    install(DIRECTORY ${CONF_DOC_DIR}/html
          DESTINATION ${CPPQED_DOC_INSTALL_DIR}/${PROJECT_NAME}
          OPTIONAL
    )
  endif()
endif()
