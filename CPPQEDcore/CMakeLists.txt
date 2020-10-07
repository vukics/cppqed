
#! \ingroup Main
#! \file
#! \brief Top level %CMake file for the C++QED core component.
#!
#! The file has the following structure:

project (cppqedcore)

message(STATUS "Configuring CPPQEDcore")

message(STATUS ${PUBLIC_LIBS})

#! \name Project options
#! @{

#! Override the git sha commit information with this value. Used in automated builds
#! where git is not available.
set(GIT_SHA_OVERRIDE "" CACHE STRING "Override the git sha version information.")

#! @}

#! \file
#! <!--#########################################################-->
#! ### Installation directories
#! <!--#########################################################-->
#!
#! This controls into which sub-directories to put %CMake config files
#! and include files when installing.

#! \name Project variables
#! @{

#! Sub-directory of `CMAKE_INSTALL_LIBDIR` into which %CMake files are installed.
set(CPPQED_CMAKE_SUBDIR "cmake/CPPQED-${CPPQED_ID}")
#! Sub-directory of `CMAKE_INSTALL_INCLUDEDIR` into which header files are installed.
set(CPPQED_INCLUDE_DIR "CPPQED-${CPPQED_ID}")
set(CPPQED_INCLUDE_SUBDIR "${CPPQED_INCLUDE_DIR}/core")
if(CPPQED_MONOLITHIC)
  set(CPPQED_INCLUDE_DIR ${CPPQED_INCLUDE_DIR} PARENT_SCOPE)
  set(CPPQED_INCLUDE_SUBDIR ${CPPQED_INCLUDE_SUBDIR} PARENT_SCOPE)
endif()

#! @}



#! \file
#! <!--#########################################################-->
#! ### Library detection
#! <!--#########################################################-->
#!
#! In this section we look for required and optional dependencies.

#! \name Project variables
#! @{

#! \brief Dependency libraries of C++QED core which should become direct dependencies of clients of core.
#!
#! Obviously this has to include all libraries which contain templates.
set(PUBLIC_LIBS)
#! \brief Dependency libraries of C++QED core which are invisible to clients of core.
#!
#! Only libraries which are linked in to core completely (e.g. GSL) can be added to this set. Template
#! libraries have to become direct dependencies of clients using core.
set(PRIVATE_LIBS)

#! @}

find_package(GSL 2.5 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

#! \name Project options
#! @{

#! Switch for boost bzip2 support
option(COMPRESSION "Boost bzip2 support" ON)

#! @}

# Boost
find_package(Boost 1.71 REQUIRED COMPONENTS iostreams)

set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})



#! \file
#! <!--#########################################################-->
#! ### Compilation of dependencies
#! <!--#########################################################-->
#!
#! This is a convenience feature to download and compile some selected
#! dependencies automatically. It uses %CMake's
#! [ExternalProject](http://www.cmake.org/cmake/help/v2.8.12/cmake.html#module:ExternalProject).
#!
#! We use a trick to install the bundled libraries. If the bundled version of a library is used,
#! we link it privately (\ref CMake::PRIVATE_LIBS instead of \ref CMake::PUBLIC_LIBS). Then
#! we add the library manually to `CPPQED_LIBRARIES` in CPPQEDConfig.cmake.in, with the correct
#! location for the build tree and installed tree, respectively.
#!
#! If this option is set, %CMake downloads and compiles the blitz++ dependency automatically.
#!
#! Takes care also of preparing the blitz (& FLENS) dependencies in case their source is distributed together with cppqed in a package ('external' folder)




# Doxygen
find_package(Doxygen 1.8)
if(DOXYGEN_FOUND AND NOT DOXYGEN_DOT_FOUND)
    set(DOXYGEN_FOUND 0)
    message(STATUS "graphviz is needed to build the documentation, disabling Doxylink")
endif()

set(PRIVATE_LIBS ${PRIVATE_LIBS} ${GSL_LIBRARIES})

#! \file
#! <!--#########################################################-->
#! ### Compiler definitions and `config.h`
#! <!--#########################################################-->
#!
#! In this section a `config.h` file is generated and saved in the build directory.
#! This header files contains preprocessor macros indicating whether FLENS and
#! boost serialization is available. C++QED source files can then conditionally
#! compile code that depends on these features. The `config.h` file only has to
#! be included where it is needed, which is an advantage over using `-D` compiler flags
#! (less code has to be recompiled).
#!
#! Also in this section, the flag `-DBOOST_RESULT_OF_USE_TR1` is set globally if needed.
#! Note that compiler warnings are set as part of CPPQED_SETUP(), which is called in the next section.

# Generate config.h
configure_file(${CPPQED_CMAKE_MODULE_PATH}/config.h.in ${PROJECT_NAME}_config.h)

#! \file
#! <!--#########################################################-->
#! ### Compilation
#! <!--#########################################################-->
#!
#! This does some initial setup (c.f. CPPQED_SETUP() and generate_version_files()) and then builds
#! the source files in the various sub-directories. This is done by trivial `CMakeLists.txt` files
#! which only contain calls to create_object_target(). This function also handles the include
#! dependencies between directories.
#!
#! Then the core library is linked and some properties like VERSION and SOVERSION are set.
#! All header files are registered for installation with the help of the function gather_includes().

CPPQED_SETUP()

generate_version_files()

include_directories(${PROJECT_BINARY_DIR}) # for generated config files

# build libC++QEDcore.so
set(CORE_SOURCE_DIRS utils quantumdata quantumoperator structure quantumtrajectory composites high)
set(quantumdata_NEEDS utils)
set(quantumoperator_NEEDS quantumdata structure utils)
set(structure_NEEDS quantumdata utils)
set(quantumtrajectory_NEEDS structure quantumdata utils)
set(composites_NEEDS structure quantumdata utils)
set(high_NEEDS composites quantumtrajectory structure quantumoperator quantumdata utils)
foreach(d ${CORE_SOURCE_DIRS})
  add_subdirectory(${d})
endforeach(d)
set(CPPQEDLIB C++QED-${CPPQED_ID})

gather_includes(CORE_SOURCE_DIRS)
set(GENERATED_SOURCE  ${CPPQED_CMAKE_MODULE_PATH}/config.h.in
                      ${CPPQED_CMAKE_MODULE_PATH}/version.h.in
                      ${CPPQED_CMAKE_MODULE_PATH}/version.cc.in
)
add_library(${CPPQEDLIB} 
    SHARED  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.cc 
            ${OBJ_TARGETS} ${core_PUBLIC_HEADERS} ${GENERATED_SOURCE})
target_link_libraries(${CPPQEDLIB} LINK_PUBLIC ${PUBLIC_LIBS} LINK_PRIVATE ${PRIVATE_LIBS})
set_target_properties(${CPPQEDLIB} PROPERTIES
      PUBLIC_HEADER "${core_PUBLIC_HEADERS}"
      INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
      VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
      SOVERSION ${CPPQED_ABI_MAJOR}
)

#! \file
#! <!--#########################################################-->
#! ### Documentation
#! <!--#########################################################-->
#!
#! Call to cppqed_documentation().

if(CPPQED_MONOLITHIC)
  cppqed_documentation(core_ "")
endif()

#! \file
#! <!--#########################################################-->
#! ### Installation
#! <!--#########################################################-->
#!
#! This section has two tasks: prepare the build tree so that it can be found by
#! other projects which have C++QED core as a dependency, and to install all required
#! files to the system.

install(TARGETS ${CPPQEDLIB}
        EXPORT CPPQEDcoreTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR}
        COMPONENT shlib
)

# Add all targets to the build-tree export set
export(TARGETS ${CPPQEDLIB}
  FILE "${PROJECT_BINARY_DIR}/CPPQEDcoreTargets.cmake")

#! \name Project options
#! @{
#
#! Enable or disable %CMake registry (c.f. \ref cmake_find_components).
option(REGISTRY "Register build trees in the cmake registry so that other projects can find them." ON)
#! @}
if(REGISTRY)
  export(PACKAGE CPPQED)
endif(REGISTRY)

# Create the CPPQEDConfig.cmake
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_BINARY_DIR}")
foreach(d ${CORE_SOURCE_DIRS})
  set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/${d}) 
endforeach(d)
set(CONF_SPECIAL_LIBRARIES ${blitz_LIBRARIES})
if(BUNDLED_BLITZ)
  set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS})
endif()
if(FLENS AND flens_FOUND)
  if(BUNDLED_FLENS)
    set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${flens_INCLUDE_DIRS})
  endif()
endif()
set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}/CPPQEDConfig.cmake"
  INSTALL_DESTINATION  "${PROJECT_BINARY_DIR}"
  PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR CONF_SPECIAL_LIBRARIES
)
write_basic_package_version_file(${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake 
  VERSION ${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}
  COMPATIBILITY ExactVersion
)
foreach(c CPPQEDUse.cmake ElementsTemplateConfig.cmake.in 
          GetGitRevisionDescription.cmake GetGitRevisionDescription.cmake.in 
          version.cc.in version.h.in component_versions.cc.in component_versions.h.in)
  configure_file(cmake/Modules/${c} ${PROJECT_BINARY_DIR}/${c} COPYONLY)
endforeach()

# ... and for the installation tree
set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_SUBDIR})
set(CONF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR})
if(BUNDLED_BLITZ OR BUNDLED_FLENS)
  set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${BUNDLED_HEADERS})
endif()
if(BUNDLED_BLITZ)
  set(CONF_SPECIAL_LIBRARIES ${BUNDLED_BLITZ_INSTALLED_LIBRARY})
else()
  set(CONF_SPECIAL_LIBRARIES ${blitz_LIBRARIES})
endif()
if(FLENS AND flens_FOUND)
  if(BUNDLED_FLENS)
    set(CONF_SPECIAL_LIBRARIES ${CONF_SPECIAL_LIBRARIES} ${BUNDLED_FLENS_INSTALLED_LIBRARY})
  else()
    set(CONF_SPECIAL_LIBRARIES ${CONF_SPECIAL_LIBRARIES} ${flens_LIBRARIES})
  endif()
endif()
configure_package_config_file(CPPQEDConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  INSTALL_DESTINATION "${CONF_CMAKE_DIR}"
  PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR CONF_SPECIAL_LIBRARIES
)

# Install the CPPQEDConfig.cmake and CPPQEDConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQEDConfig.cmake"
  "${PROJECT_BINARY_DIR}/CPPQEDConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}/CPPQEDUse.cmake"
  "${PROJECT_BINARY_DIR}/ElementsTemplateConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/GetGitRevisionDescription.cmake"
  "${PROJECT_BINARY_DIR}/GetGitRevisionDescription.cmake.in"
  "${PROJECT_BINARY_DIR}/version.h.in"
  "${PROJECT_BINARY_DIR}/version.cc.in"
  "${PROJECT_BINARY_DIR}/component_versions.h.in"
  "${PROJECT_BINARY_DIR}/component_versions.cc.in"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT CPPQEDcoreTargets DESTINATION
  "${CMAKE_INSTALL_LIBDIR}/${CPPQED_CMAKE_SUBDIR}" COMPONENT dev)


#! \file
#! <!--#########################################################-->
#! ### Installation of dependencies
#! <!--#########################################################-->
#!
#! If blitz and flens are compiled as external projects, install them as well. Libraries
#! and header files are installed into the same locations as CPPQEDcore, this way it is
#! ensured they are found.

set(BUNDLED_HEADERS "${CMAKE_INSTALL_INCLUDEDIR}/${CPPQED_INCLUDE_DIR}/bundled")

if(BUNDLED_BLITZ)
  ExternalProject_Get_Property(blitz install_dir)
  file(GLOB blitz_libs ${install_dir}/lib/*.so*)
  set(blitz_libs ${install_dir}/lib/libblitz.so
                 ${install_dir}/lib/libblitz.so.0
                 ${install_dir}/lib/libblitz.so.1.0.3) # it is a bit inconvenient to have to write the version explicitly, but fortunately it doesn’t change very often
  install(FILES ${blitz_libs}
          PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${install_dir}/include/blitz
          DESTINATION ${BUNDLED_HEADERS})
endif()

if(BUNDLED_FLENS)
  ExternalProject_Get_Property(flens install_dir)
  install(DIRECTORY ${install_dir}/include/flens
          DESTINATION ${BUNDLED_HEADERS})
endif()


#! \file
#! <!--#########################################################-->
#! ### Compile extras
#! <!--#########################################################-->
#!
#! Compile the example code documented \ref structurebundleguide "here".

add_subdirectory(examples)

#! \file
#! <!--#########################################################-->
#! ### Summary of enabled/disabled features
#! <!--#########################################################-->
#!
#! Display some nice summary of which components have been found and which not.

set_package_properties(PkgConfig PROPERTIES URL "http://pkgconfig.freedesktop.org/wiki"
                                DESCRIPTION "Package config system that manages compile/link flags"
                                TYPE OPTIONAL
                                PURPOSE "Assists cmake in finding libraries.")
set_package_properties(blitz PROPERTIES  URL "http://github.com/blitzpp/blitz"
                                DESCRIPTION "High-performance C++ vector mathematics library"
                                TYPE REQUIRED
                                PURPOSE "Multi-Array implementation used in the framework.")
set_package_properties(flens PROPERTIES URL "http://flens.sourceforge.net"
                                DESCRIPTION "Flexible Library for Efficient Numerical Solutions."
                                TYPE OPTIONAL
                                PURPOSE "With FLENS the negativity of the partially transposed density operator can be calculated.")
set_package_properties(Boost PROPERTIES URL "http://www.boost.org/"
                                DESCRIPTION "Collection of portable C++ source libraries."
                                TYPE REQUIRED
                                PURPOSE "Advanced template metaprogramming and preprocessor algorithms.")
set_package_properties(GSL PROPERTIES URL "http://www.gnu.org/software/gsl/"
                                DESCRIPTION "GNU Scientific Library"
                                TYPE REQUIRED
                                PURPOSE "Used as implementation of an ODE solver and random number generator.")
set_package_properties(Doxygen PROPERTIES URL "http://www.doxygen.org/"
                                DESCRIPTION "Generate documentation from source code."
                                TYPE OPTIONAL
                                PURPOSE "Generation of API documentation. Note that you also need the graphviz package to build the documentation.")
add_feature_info("FLENS" CPPQED_FLENS_FOUND "compile framework with FLENS support.")
add_feature_info("Compression" COMPRESSION "bzip2 compression of binary statevectors.")
if(NOT ${CPPQED_MONOLITHIC})
  feature_summary( WHAT ALL )
endif()

