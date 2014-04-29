
#! \ingroup Main
#! \file
#! \brief Top level %CMake file for the C++QED core component.
#!
#! The file has the following structure:

cmake_minimum_required (VERSION 2.8.9)

project (core)

#! \name Project variables
#! @{

#! \brief Path to additional CMake modules.
set(CPPQED_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CPPQED_CMAKE_MODULE_PATH})

#! @}

include(FeatureSummary)
include(CPPQEDUse)

#! \file
#! <!--#########################################################-->
#! ### Version management
#! <!--#########################################################-->
#!
#! Set the major, minor and patch version of C++QED (c.f. \ref versioning). The major and minor
#! versions end up in the library names of all components. In this section,
#! also the ABI version of the library is set (c.f. \ref ABI).


#! \anchor versioning
#! \name Version variables
#! Versioning scheme
#! @{

#! MAJOR - the grand version of C++QED, as in v2
set(CPPQED_VERSION_MAJOR 2)
#! MINOR - the milestone release, as in v2m10
set(CPPQED_VERSION_MINOR 100)
#! PATCH - should be increased on every snapshot package release, reset on milestone change
set(CPPQED_VERSION_PATCH 2)
#! Full C++QED version
set(CPPQED_VERSION
  "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}.${CPPQED_VERSION_PATCH}")
#! MAJOR.MINOR, at the moment mainly distinguishes stable (2.10) from development (2.100)
set(CPPQED_ID "${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}")

#! @}

#! \anchor ABI
#! \name ABI versioning scheme
#!
#! Adopt the shared library versioning scheme of libtool
#! - CURRENT  - the current ABI version
#! - AGE      - number of versions backwards that CURRENT is compatible with
#! - REVISION - patch revision
#!
#! *Rationale*
#!
#! Bugfixes and Patches which do not affect ABI:
#!  * increase REVISION
#!
#! Symbols added to library (i.e. binary compatibility NOT broken):
#!  * increase AGE
#!  * remember to call dh_makeshlibs with -V packagename (>=packageversion) when packaging
#!
#! Symbols removed or changed (i.e. binary compatibility broken):
#!  * set CURRENT -> CURRENT + AGE + 1, reset AGE
#!
#! On every change of AGE or CURRENT: reset REVISION
#!
#! The library version (.so.X.Y.Z) is `{CURRENT-AGE}.{AGE}.{REVISION}`, here named
#! `{ABI_MAJOR}.{ABI_MINOR}.{ABI_MICRO}`. This way it is ensured that ABI_MAJOR only changes
#! when binary compatibility is broken.
#!
#! The SONAME of the library always is: `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}.so.${ABI_MAJOR}`
#! and the packages are named `libC++QED-${MAJOR_VERSION}.${MINOR_VERSION}-${ABI_MAJOR}[-dev]`
#! @{

#! The current ABI version
set(CPPQED_ABI_CURRENT 3)
#! number of Versions backwards that CURRENT is compatible with
set(CPPQED_ABI_AGE 0)
#! Patch revision
set(CPPQED_ABI_REVISION 0)
#! Set to ABI_CURRENT - ABI_AGE
set(CPPQED_ABI_MAJOR)
math(EXPR CPPQED_ABI_MAJOR "${CPPQED_ABI_CURRENT}-${CPPQED_ABI_AGE}")
#! Same as ABI_AGE
set(CPPQED_ABI_MINOR "${CPPQED_ABI_AGE}")
#! Same as ABI_REVISION
set(CPPQED_ABI_MICRO "${CPPQED_ABI_REVISION}")

#! @}

if(CPPQED_MONOLITHIC)
  set(CPPQED_VERSION ${CPPQED_VERSION} PARENT_SCOPE)
  set(CPPQED_VERSION_MAJOR ${CPPQED_VERSION_MAJOR} PARENT_SCOPE)
  set(CPPQED_VERSION_MINOR ${CPPQED_VERSION_MINOR} PARENT_SCOPE)
  set(CPPQED_VERSION_PATCH ${CPPQED_VERSION_PATCH} PARENT_SCOPE)
  set(CPPQED_ID ${CPPQED_ID} PARENT_SCOPE)
  set(CPPQED_ABI_MAJOR ${CPPQED_ABI_MAJOR} PARENT_SCOPE)
endif()

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
#! ### Compiler detection
#! <!--#########################################################-->
#!
#! At the moment g++ >= 4.7 and Clang >= 3.1 is needed for C++11 features.

set(GCC_MINIMAL 4.7)
set(CLANG_MINIMAL 3.1)


if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${GCC_MINIMAL})
  message(FATAL_ERROR "GCC g++ version >= ${GCC_MINIMAL} needed.")
endif ()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MINIMAL})
  message(FATAL_ERROR "Clang version >= ${CLANG_MINIMAL} needed.")
endif ()

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

find_package(GSL REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

#! \name Project options
#! @{

#! Switch for boost serialization.
option(SERIALIZATION "Boost serialization" ON)
if(SERIALIZATION)
  set(CPPQED_SERIALIZATION_FOUND 1)
else(SERIALIZATION)
  set(CPPQED_SERIALIZATION_FOUND 0)
endif(SERIALIZATION)

#! Switch for FLENS support.
option(FLENS "FLENS support" ON)

#! @}

# blitz
if(NOT BUNDLED_BLITZ)
  find_package(blitz REQUIRED)
  set(PUBLIC_LIBS ${PUBLIC_LIBS} ${blitz_LIBRARIES})
  set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS})
else()
  message(STATUS "Using bundled blitz from ${blitz_LIBRARIES}.")
  set(PRIVATE_LIBS ${PRIVATE_LIBS} ${blitz_LIBRARIES})
  include_directories(SYSTEM ${blitz_INCLUDE_DIRS})
endif()
if( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built with serialization support.")
else( blitz_SERIALIZATION_FOUND )
  message(STATUS "Blitz++ built without serialization support. Please configure blitz with --enable-serialization to enable optional binary statevector output.")
  set(CPPQED_SERIALIZATION_FOUND 0)
endif( blitz_SERIALIZATION_FOUND )

# Boost (OPTIONAL_COMPONENTS does not work with Boost find_package)
find_package(Boost 1.46.0 REQUIRED)
find_package(Boost QUIET COMPONENTS serialization)

if( Boost_SERIALIZATION_FOUND )
  message(STATUS "Boost serialization library found.")
else( Boost_SERIALIZATION_FOUND )
  message(STATUS "Boost serialization library not found.")
  set(CPPQED_SERIALIZATION_FOUND 0)
endif( Boost_SERIALIZATION_FOUND )
set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})


# FLENS
if(NOT BUNDLED_FLENS)
  find_package(flens)
endif()
if( flens_FOUND AND FLENS )
  include_directories(SYSTEM ${flens_INCLUDE_DIRS})
  if(BUNDLED_FLENS)
    include_directories(SYSTEM ${flens_INCLUDE_DIRS})
    message(STATUS "Using bundled flens from ${flens_LIBRARIES}.")
  else()
    set(CPPQED_THIRDPARTY_INCLUDE_DIRS ${CPPQED_THIRDPARTY_INCLUDE_DIRS} ${flens_INCLUDE_DIRS})
  endif()
  set(CPPQED_FLENS_FOUND 1)
else( flens_FOUND AND FLENS )
  message(STATUS "Flens library not found or disabled, optional flens support disabled.")
  set(DO_NOT_USE_FLENS ON)
  set(CPPQED_FLENS_FOUND 0)
endif( flens_FOUND AND FLENS )

# Check if serialization can be enabled, inform user
if( CPPQED_SERIALIZATION_FOUND )
  message(STATUS "Support for binary statevector output enabled.")
  set(PUBLIC_LIBS ${PUBLIC_LIBS} ${Boost_SERIALIZATION_LIBRARY})
else( CPPQED_SERIALIZATION_FOUND )
  set(DO_NOT_USE_BOOST_SERIALIZATION ON)
  message(STATUS "Optional support for binary statevector output disabled.")
endif( CPPQED_SERIALIZATION_FOUND )

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

if(Boost_FOUND)
  if(Boost_MAJOR_VERSION EQUAL "1" AND Boost_MINOR_VERSION GREATER "51")
    set(CPPQED_DEFINITIONS ${CPPQED_DEFINITIONS} -DBOOST_RESULT_OF_USE_TR1)
    message(STATUS "added -DBOOST_RESULT_OF_USE_TR1" )
  endif()
endif()


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
set(CORE_SOURCE_DIRS utils quantumdata quantumoperator structure quantumtrajectory composites)
set(quantumdata_NEEDS utils)
set(quantumoperator_NEEDS quantumdata structure utils)
set(structure_NEEDS quantumdata utils)
set(quantumtrajectory_NEEDS structure quantumdata utils)
set(composites_NEEDS structure quantumdata utils)
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
  set(CONF_SPECIAL_LIBRARIES ${CONF_SPECIAL_LIBRARIES} ${flens_LIBRARIES})
endif()
set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
set(CONF_FLAT_INCLUDE 0)
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
set(CONF_FLAT_INCLUDE 1)
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
set_package_properties(blitz PROPERTIES  URL "http://sf.net/projects/blitz/"
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
                                DESCRIPTION "Generate documentation from source code"
                                TYPE OPTIONAL
                                PURPOSE "Generation of API documentation.")
add_feature_info("FLENS" CPPQED_FLENS_FOUND "compile framework with FLENS support.")
add_feature_info("Serialization" CPPQED_SERIALIZATION_FOUND "needed for binary statevector output." )
if(NOT ${CPPQED_MONOLITHIC})
  feature_summary( WHAT ALL )
endif()

