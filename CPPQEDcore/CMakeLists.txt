# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring CPPQEDcore")

add_subdirectory(quantumdata)

add_subdirectory(structure)

add_subdirectory(quantumoperator)

add_subdirectory(quantumtrajectory)

add_subdirectory(composites)

add_subdirectory(high)

add_library(CPPQEDcore SHARED
            $<TARGET_OBJECTS:CPPQED::CPPQEDcore_quantumdata> $<TARGET_OBJECTS:CPPQED::CPPQEDcore_structure>
            $<TARGET_OBJECTS:CPPQED::CPPQEDcore_quantumoperator> $<TARGET_OBJECTS:CPPQED::CPPQEDcore_quantumtrajectory>
            $<TARGET_OBJECTS:CPPQED::CPPQEDcore_composites> $<TARGET_OBJECTS:CPPQED::CPPQEDcore_high>)
add_library(CPPQED::CPPQEDcore ALIAS CPPQEDcore)

target_link_libraries(CPPQEDcore PUBLIC CPPQED::CPPQEDcore_high)

install(TARGETS CPPQEDcore DESTINATION lib)

## Doxygen
#find_package(Doxygen 1.8)
#if(DOXYGEN_FOUND AND NOT DOXYGEN_DOT_FOUND)
    #set(DOXYGEN_FOUND 0)
    #message(STATUS "graphviz is needed to build the documentation, disabling Doxylink")
#endif()


## build libC++QEDcore.so
#set(CORE_SOURCE_DIRS utils quantumdata quantumoperator structure quantumtrajectory composites high)
#set(quantumdata_NEEDS utils)
#set(structure_NEEDS quantumdata utils)
#set(quantumoperator_NEEDS quantumdata structure utils)
#set(quantumtrajectory_NEEDS structure quantumdata utils)
#set(composites_NEEDS structure quantumdata utils)
#set(high_NEEDS composites quantumtrajectory structure quantumoperator quantumdata utils)
#foreach(d ${CORE_SOURCE_DIRS})
  #add_subdirectory(${d})
#endforeach(d)
#set(CPPQEDLIB C++QED-${CPPQED_ID})

#gather_includes(CORE_SOURCE_DIRS)
#set(GENERATED_SOURCE  ${CPPQED_CMAKE_MODULE_PATH}/config.h.in
                      #${CPPQED_CMAKE_MODULE_PATH}/version.h.in
                      #${CPPQED_CMAKE_MODULE_PATH}/version.cc.in
#)
#add_library(${CPPQEDLIB} 
    #SHARED  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.cc 
            #${OBJ_TARGETS} ${core_PUBLIC_HEADERS} ${GENERATED_SOURCE})
#target_link_libraries(${CPPQEDLIB} LINK_PUBLIC ${PUBLIC_LIBS} LINK_PRIVATE ${PRIVATE_LIBS})
#set_target_properties(${CPPQEDLIB} PROPERTIES
      #PUBLIC_HEADER "${core_PUBLIC_HEADERS}"
      #INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
      #VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
      #SOVERSION ${CPPQED_ABI_MAJOR}
#)

##! \file
##! <!--#########################################################-->
##! ### Documentation
##! <!--#########################################################-->
##!
##! Call to cppqed_documentation().

#if(CPPQED_MONOLITHIC)
  #cppqed_documentation(core_ "")
#endif()

##! \file
##! <!--#########################################################-->
##! ### Installation
##! <!--#########################################################-->
##!
##! This section has two tasks: prepare the build tree so that it can be found by
##! other projects which have C++QED core as a dependency, and to install all required
##! files to the system.

#install(TARGETS ${CPPQEDLIB}
        #EXPORT CPPQEDcoreTargets
        #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        #PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR}
        #COMPONENT shlib
#)

## Add all targets to the build-tree export set
#export(TARGETS ${CPPQEDLIB}
  #FILE "${PROJECT_BINARY_DIR}/CPPQEDcoreTargets.cmake")

##! \name Project options
##! @{
##
##! Enable or disable %CMake registry (c.f. \ref cmake_find_components).
#option(REGISTRY "Register build trees in the cmake registry so that other projects can find them." ON)
##! @}
#if(REGISTRY)
  #export(PACKAGE CPPQED)
#endif(REGISTRY)

## Create the CPPQEDConfig.cmake
## ... for the build tree
#set(CONF_INCLUDE_DIRS "${PROJECT_BINARY_DIR}")
#foreach(d ${CORE_SOURCE_DIRS})
  #set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/${d}) 
#endforeach(d)
#set(CONF_SPECIAL_LIBRARIES ${blitz_LIBRARIES})
#if(BUNDLED_BLITZ)
  #set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS})
#endif()
#set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
#configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}/CPPQEDConfig.cmake"
  #INSTALL_DESTINATION  "${PROJECT_BINARY_DIR}"
  #PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR CONF_SPECIAL_LIBRARIES
#)
#write_basic_package_version_file(${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake 
  #VERSION ${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}
  #COMPATIBILITY ExactVersion
#)
#foreach(c CPPQEDUse.cmake ElementsTemplateConfig.cmake.in 
          #GetGitRevisionDescription.cmake GetGitRevisionDescription.cmake.in 
          #version.cc.in version.h.in component_versions.cc.in component_versions.h.in)
  #configure_file(cmake/Modules/${c} ${PROJECT_BINARY_DIR}/${c} COPYONLY)
#endforeach()

## ... and for the installation tree
#set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR})
#set(CONF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR})
#if(BUNDLED_BLITZ)
  #set(CONF_SPECIAL_LIBRARIES ${BUNDLED_BLITZ_INSTALLED_LIBRARY})
#else()
  #set(CONF_SPECIAL_LIBRARIES ${blitz_LIBRARIES})
#endif()
#configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  #INSTALL_DESTINATION "${CONF_CMAKE_DIR}"
  #PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR CONF_SPECIAL_LIBRARIES
#)

## Install the CPPQEDConfig.cmake and CPPQEDConfigVersion.cmake
#install(FILES
  #"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  #"${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake"
  #"${PROJECT_BINARY_DIR}/CPPQEDUse.cmake"
  #"${PROJECT_BINARY_DIR}/ElementsTemplateConfig.cmake.in"
  #"${PROJECT_BINARY_DIR}/GetGitRevisionDescription.cmake"
  #"${PROJECT_BINARY_DIR}/GetGitRevisionDescription.cmake.in"
  #"${PROJECT_BINARY_DIR}/version.h.in"
  #"${PROJECT_BINARY_DIR}/version.cc.in"
  #"${PROJECT_BINARY_DIR}/component_versions.h.in"
  #"${PROJECT_BINARY_DIR}/component_versions.cc.in"
  #DESTINATION "${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)
 
## Install the export set for use with the install-tree
#install(EXPORT CPPQEDcoreTargets DESTINATION
  #"${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)


##! \file
##! <!--#########################################################-->
##! ### Installation of dependencies
##! <!--#########################################################-->
##!


##! \file
##! <!--#########################################################-->
##! ### Compile extras
##! <!--#########################################################-->
##!
##! Compile the example code documented \ref structurebundleguide "here".

#add_subdirectory(examples)

##! \file
##! <!--#########################################################-->
##! ### Summary of enabled/disabled features
##! <!--#########################################################-->
##!
##! Display some nice summary of which components have been found and which not.

#set_package_properties(PkgConfig PROPERTIES URL "http://pkgconfig.freedesktop.org/wiki"
                                #DESCRIPTION "Package config system that manages compile/link flags"
                                #TYPE OPTIONAL
                                #PURPOSE "Assists cmake in finding libraries.")
#set_package_properties(blitz PROPERTIES  URL "http://github.com/blitzpp/blitz"
                                #DESCRIPTION "High-performance C++ vector mathematics library"
                                #TYPE REQUIRED
                                #PURPOSE "Multi-Array implementation used in the framework.")
#set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                #DESCRIPTION "Collection of portable C++ source libraries."
                                #TYPE REQUIRED
                                #PURPOSE "Advanced template metaprogramming and preprocessor algorithms.")
#set_package_properties(GSL PROPERTIES URL "http://www.gnu.org/software/gsl/"
                                #DESCRIPTION "GNU Scientific Library"
                                #TYPE REQUIRED
                                #PURPOSE "Used as implementation of an ODE solver and random number generator.")
#set_package_properties(Doxygen PROPERTIES URL "http://www.doxygen.org/"
                                #DESCRIPTION "Generate documentation from source code."
                                #TYPE OPTIONAL
                                #PURPOSE "Generation of API documentation. Note that you also need the graphviz package to build the documentation.")
#add_feature_info("Compression" COMPRESSION "bzip2 compression of binary statevectors.")
#if(NOT ${CPPQED_MONOLITHIC})
  #feature_summary( WHAT ALL )
#endif()

