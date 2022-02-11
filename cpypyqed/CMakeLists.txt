# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)


message(STATUS "Configuring cpypyqed")

find_package(Boost ${BOOST_MINIMAL} COMPONENTS python3 numpy3)

find_package(Python 3.8 COMPONENTS Development)

if ( ( NOT TARGET Boost::python3 ) OR ( NOT TARGET Boost::numpy3 ) )
  return()
endif()


# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

add_library(cpypyqed MODULE io.cc) # MODULE means it cannot be linked to, only loaded @ runtime

add_library(CPPQED::cpypyqed ALIAS cpypyqed)

target_link_libraries(cpypyqed PRIVATE CPPQED::CPPQEDcore Boost::python3 Boost::numpy3 Python::Python)


#! \file
#! <!--#########################################################-->
#! ### Documentation
#! <!--#########################################################-->
#!
#! Find Doxygen, Sphinx and Doxylink. If all these components are found,
#! build the documentation. Sphinx has the advantage that it actually loads
#! the modules with a python interpreter and looks at the docstrings.
#! Otherwise it would be very cumbersome to document a boost python module,
#! because all the signatures would have to be duplicated in the documentation
#! code.
#!
#! The relative path to the Doxygen tagfile `core.tag` is established, in order
#! to link to the C++QED API with Doxylink.
#!
#! In monolithic mode, the file `doc`layout.html.in` is copied over to the
#! documentation build directory. This files contains links to the other C++QED
#! components and makes them appear on the Sphinx pages.
#!
#! The sphinx source file `reference.rst` depends on the Python package name,
#! therefore we again have to generate this file at build time, when the
#! build configuration is known in any case. The solution is again, as in the
#! previous section, to let a custom command generate the file with a simple
#! cmake script.
#!
#! Finally, the call to the sphinx binary is added as a custom command, the
#! target is named `cpypyqed_doc`. In monolithic mode, the resulting html
#! directory is installed to the system.

#find_package(Sphinx)
#find_package(Doxylink)

#if(CPPQED_MONOLITHIC AND DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND AND SPHINX_FOUND AND DOXYLINK_FOUND)
  #set(CPPQED_DOC_DIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/cppqed-doc-${CPPQED_ID}")
  #file(MAKE_DIRECTORY ${${PROJECT_NAME}_DOC_DIR})
  #file(RELATIVE_PATH core_RELATIVE_DOXYGEN_TAG ${${PROJECT_NAME}_DOC_DIR} ${core_DOXYGEN_TAG})
  #get_filename_component(core_RELATIVE_DOXYGEN_DIR ${core_RELATIVE_DOXYGEN_TAG} PATH)
  #configure_file(${PROJECT_SOURCE_DIR}/doc/conf.py.in ${${PROJECT_NAME}_DOC_DIR}/conf.py @ONLY)
  #if(CPPQED_MONOLITHIC)
    #configure_file(${PROJECT_SOURCE_DIR}/doc/layout.html.in ${${PROJECT_NAME}_DOC_DIR}/_templates/layout.html @ONLY)
  #endif()
  #foreach(AUXFILE index.rst userguide.rst tools.rst)
    #configure_file(${PROJECT_SOURCE_DIR}/doc/${AUXFILE} ${${PROJECT_NAME}_DOC_DIR} COPYONLY)
  #endforeach()
  ## Sphinx cache with pickled ReST documents
  #set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
  ## HTML output directory
  #set(SPHINX_HTML_DIR "${${PROJECT_NAME}_DOC_DIR}/html")
  #add_custom_command(OUTPUT ${${PROJECT_NAME}_DOC_DIR}/reference.rst
                     #COMMAND ${CMAKE_COMMAND}
                        #-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/doc/reference.rst.in
                        #-DOUTPUT=${${PROJECT_NAME}_DOC_DIR}/reference.rst
                        #-DDEBUG_SUFFIX="$<$<CONFIG:Debug>:_d>"
                        #-P ${CMAKE_CURRENT_SOURCE_DIR}/generate.cmake
                     #DEPENDS generate.cmake doc/reference.rst.in
  #)
  #add_custom_target(${PROJECT_NAME}_doc
      #${SPHINX_EXECUTABLE}
          #-q -b html
          #-c "${${PROJECT_NAME}_DOC_DIR}"
          #-d "${SPHINX_CACHE_DIR}"
          #"${${PROJECT_NAME}_DOC_DIR}"
          #"${SPHINX_HTML_DIR}"
      #WORKING_DIRECTORY ${${PROJECT_NAME}_DOC_DIR}
      #COMMENT "Building HTML documentation with Sphinx"
      #DEPENDS ${${PROJECT_NAME}_DOC_DIR}/reference.rst
  #)
  #add_dependencies(${PROJECT_NAME}_doc cpypyqed)
  #if(CPPQED_MONOLITHIC)
    #add_dependencies(${PROJECT_NAME}_doc core_doc elements_doc)
    #install(DIRECTORY ${${PROJECT_NAME}_DOC_DIR}/html
          #DESTINATION ${CPPQED_DOC_DIR}/${PROJECT_NAME}
          #OPTIONAL
    #)
  #endif()
#endif()

#! \file
#! <!--#########################################################-->
#! ### Feature summary
#! <!--#########################################################-->
#!
#! Generate a summary of the cpypyqed dependencies.

#set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                #DESCRIPTION "Collection of portable C++ source libraries."
                                #TYPE REQUIRED
                                #PURPOSE "Python wrapper library for C++ libraries.")
#set_package_properties(PythonLibs PROPERTIES URL "http://www.python.org/"
                                #DESCRIPTION "General-purpose, high-level programming language."
                                #TYPE OPTIONAL
                                #PURPOSE "Access C++QED functionality through Python (work in progress, incomplete).")
#set_package_properties(Numpy PROPERTIES URL "http://www.numpy.org/"
                                #DESCRIPTION "NumPy is the fundamental package for scientific computing with Python"
                                #TYPE OPTIONAL
                                #PURPOSE "Data input/output with Python, testsuite.")
#set_package_properties(Scipy PROPERTIES URL "http://www.scipy.org/"
                                #DESCRIPTION "SciPy is a Python-based ecosystem of open-source software for mathematics, science, and engineering."
                                #TYPE OPTIONAL
                                #PURPOSE "The testsuite depends on this.")
#set_package_properties(Sphinx PROPERTIES URL "http://sphinx-doc.org/"
                                #DESCRIPTION "Sphinx is a documentation generator which converts reStructuredText files into HTML websites"
                                #TYPE OPTIONAL
                                #PURPOSE "Sphinx is needed to build the cpypyqed documentation.")
#set_package_properties(Doxylink PROPERTIES URL "https://pypi.python.org/pypi/sphinxcontrib-doxylink"
                                #DESCRIPTION "Doxylink is a Sphinx extension to link to external Doxygen API documentation."
                                #PURPOSE "Needed to link to the rest of the C++QED in the cpypyqed documentation.")

#if(NOT CPPQED_MONOLITHIC)
  #feature_summary( WHAT ALL )
#endif()
