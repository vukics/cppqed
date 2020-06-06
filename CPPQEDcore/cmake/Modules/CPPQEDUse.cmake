# Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \file CPPQEDUse.cmake
#! \brief Macros and functions which help to build C++QED projects.

#! \addtogroup CPPQEDUse
#!  @{

#! \brief Set common compiler flags for C++QED projects.
#!
#! \return This function adds flags to the `CMAKE_CXX_FLAGS` and the
#! `CMAKE_CXX_FLAGS_DEBUG` variables.
#!
#! This macro sets `-std=c++11` and otherwise only affects warnings. In Debug mode, warnings are enabled.
#! It is automatically called in CPPQED_SETUP().
macro(cppqed_cxx_flags)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-ignored-qualifiers -Wno-sign-compare -Wno-overloaded-virtual")
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-local-type-template-args")
   if(CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3)
     if(XCODE)
       set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
     else()
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --stdlib=libc++")
     endif()
   endif()
  endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
endmacro(cppqed_cxx_flags)

#! \brief Generate a list of all header files in the current project.
#! \param source_dirs <var-name> Variable name of a list containing all source directories (e.g. `SOURCE_DIRS`, not `${SOURCE_DIRS}`)
#! \return This function sets the variable `${PROJECT_NAME}_PUBLIC_HEADERS`.
#!
#! Typically the output of this function is used in the `PUBLIC_HEADER` property of libraries and install targets.
function(gather_includes source_dirs)
  foreach(d ${${source_dirs}})
    file(GLOB INC ${d}/*.h ${d}/*.tcc)
    set(${PROJECT_NAME}_PUBLIC_HEADERS ${${PROJECT_NAME}_PUBLIC_HEADERS} ${INC})
  endforeach(d)
  foreach(GENERATED_HEADER ${PROJECT_NAME}_config.h ${PROJECT_NAME}_version.h component_versions.h)
    if(EXISTS ${PROJECT_BINARY_DIR}/${GENERATED_HEADER})
      set(${PROJECT_NAME}_PUBLIC_HEADERS ${${PROJECT_NAME}_PUBLIC_HEADERS} ${PROJECT_BINARY_DIR}/${GENERATED_HEADER})
    endif()
  endforeach()
  set(${PROJECT_NAME}_PUBLIC_HEADERS ${${PROJECT_NAME}_PUBLIC_HEADERS} PARENT_SCOPE)
endfunction()


#! \brief Gather all source files in the current directory and create an object target.
#! \return This macro sets the variable `OBJ_TARGETS`, which can be used as
#!    source for library targets.
#!
#! This macro creates a target `${PROJECT_NAME}_${name}_objs` where `${name}` is the current
#! subdirectory and populates it with all source files in this directory. It then adds this
#! target to the parent scope variable `OBJ_TARGETS`.
#!
#! The macro processes the list `${name}_NEEDS`. This list contains the name of all the other
#! subdirectories (relative to the top level) from which this subdirectory may include header files.
macro(create_object_target)
  get_filename_component(name ${CMAKE_CURRENT_LIST_DIR} NAME)
  include_directories(${CMAKE_CURRENT_LIST_DIR})
  foreach(d ${${name}_NEEDS})
    include_directories(${PROJECT_SOURCE_DIR}/${d})
  endforeach(d)
  aux_source_directory(. ${name}_srcs)
  add_library(${PROJECT_NAME}_${name}_objs OBJECT ${${name}_srcs})
  if(BUNDLED_BLITZ)
    add_dependencies(${PROJECT_NAME}_${name}_objs cppqed-blitz)
  endif()
  if(BUNDLED_FLENS)
    add_dependencies(${PROJECT_NAME}_${name}_objs flens)
  endif()
  set_target_properties(${PROJECT_NAME}_${name}_objs PROPERTIES POSITION_INDEPENDENT_CODE On)
  set(OBJ_TARGETS ${OBJ_TARGETS} "\$<TARGET_OBJECTS:${PROJECT_NAME}_${name}_objs>" PARENT_SCOPE)
endmacro()

if(NOT ${CPPQED_MONOLITHIC})
  message(STATUS "Using CPPQED_SETUP macro from ${CMAKE_CURRENT_LIST_DIR}")
endif()

#! \brief Initialize a C++QED project.
#!
#! Typical usage:
#!
#!     find_package(CPPQED 2.100 REQUIRED)
#!     include(${CPPQED_USE})
#!     CPPQED_SETUP()
#!
#! This macro makes sure all C++QED sub-projects are initialized in the same way.
#! Installation directories are handled by `GNUInstallDirs`, which is automatically included here.
#! The following actions are performed:
#!
#! - Sets the build type to Release if none was specified with `-DCMAKE_BUILD_TYPE`. If an unknown
#!  build type is encountered this will result in an error.
#! - Sets appropriate compiler flags by calling cppqed_cxx_flags() and adding `-DBZ_DEBUG` as
#!  flag in debug mode, adding ::CPPQED_DEFINITIONS to compiler flags.
#! - Adds C++QED core include directories and third-party include directories (boost, blitz etc.)
#!  to the include path.
#! - Handles the libraries rpath: When building, add the full rpath into the libraries which is needed
#!  to find all dependencies (also outside the build tree). When installing the libraries, keep only those
#!  paths which are not system paths. See [here](http://www.cmake.org/Wiki/CMake_RPATH_handling) for more
#!  information.
macro(CPPQED_SETUP)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)
  if(DEFINED CPPQED_CMAKE_DIR)
    include(${CPPQED_CMAKE_DIR}/GetGitRevisionDescription.cmake)
  else()
    include(GetGitRevisionDescription)
  endif()

  # guard against bad build-type strings

  if (NOT CMAKE_BUILD_TYPE)
    message(WARNING "Build type not set, default is \"Release\".")
    set(CMAKE_BUILD_TYPE "Release")
  endif()

  string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
  if(   NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug and Release (case-insensitive).")
  endif()

  set(CMAKE_DEBUG_POSTFIX "_d")

  cppqed_cxx_flags()

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG")
  add_definitions(${CPPQED_DEFINITIONS})

  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${CMAKE_INSTALL_FULL_LIBDIR} isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})
  ENDIF("${isSystemDir}" STREQUAL "-1")

  include_directories(${CPPQED_INCLUDE_DIRS})
  include_directories(SYSTEM ${CPPQED_THIRDPARTY_INCLUDE_DIRS})

endmacro()

#! \brief Initialize an elements project.
#!
#! This macro expects `ELEMENTS_SOURCE_DIRS` to be set to a list of source directories. Include dependencies between
#! subdirectories can be established by setting `<subdir>_NEEDS` variables.
#!
#! It is used for both the original C++QED elements project and custom element projects. The name of a custom
#! elements project may not be "elements", this is checked here. The following actions are performed:
#!
#! - Find C++QED core, and additionally for custom element projects find C++QED elements and add them as dependencies.
#! - Generate version information to compile into the libraries, see generate_version_files() for details.
#! - Process all subdirectories listed in `ELEMENTS_SOURCE_DIRS`. The `CMakeLists.txt` file in the subdirectories should
#!  only have a single call to create_object_target().
#! - Create the library target `C++QED${PROJECT_NAME}-${CPPQED_ID}`,  e.g. `C++QEDelements-2.100` and link to all dependencies.
#! - The version and SONAME is the same as for the C++QED core library.
#! - Create an install target to install the libraries and the headers.
#! - Generate an appropriate `CPPQED${PROJECT_NAME}ConfigVersion.cmake` and `CPPQED${PROJECT_NAME}Config.cmake` file, which
#!  ensure that other projects can find the current library by calling, e.g. `find_package(CPPQEDelements)`. Two versions of these
#!  files are generated, one suitable for the build tree and one which is getting installed along with the library.
#!  For details, see \ref Config "this" page.
#!
#! Typical usage in a custom elements project:
#!
#!     project(elements_custom)
#!     find_package(CPPQED 2.100 REQUIRED)
#!     include(${CPPQED_USE})
#!
#!     set(ELEMENTS_SOURCE_DIRS utils frees interactions)
#!     set(frees_NEEDS utils)
#!     set(interactions_NEEDS utils frees)
#!
#!     elements_project()
macro(elements_project)
  if(${PROJECT_NAME} STREQUAL elements)
    if(NOT DEFINED ORIGINAL_ELEMENTS_PROJECT)
      message(FATAL_ERROR "The project cannot be named 'elements', as this is the name of the CPPQED elements project.")
    endif()
  else()
    find_package(CPPQEDelements ${CPPQED_ID} REQUIRED)
    include_directories(${CPPQEDelements_INCLUDE_DIRS})
    message(STATUS "Using C++QED elements from ${CPPQEDelements_DIR}.")
  endif()

  CPPQED_SETUP()

  generate_version_files()

  foreach(d ${ELEMENTS_SOURCE_DIRS})
    add_subdirectory(${d})
  endforeach(d)

  set(ELEMENTS_CMAKE_SUBDIR "cmake/CPPQED${PROJECT_NAME}-${CPPQED_ID}")
  set(ELEMENTS_INCLUDE_SUBDIR "CPPQED-${CPPQED_ID}/${PROJECT_NAME}")

  gather_includes(ELEMENTS_SOURCE_DIRS)

  set(ELEMENTSLIB C++QED${PROJECT_NAME}-${CPPQED_ID})
  add_library(${ELEMENTSLIB} SHARED ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.cc ${${PROJECT_NAME}_PUBLIC_HEADERS} ${OBJ_TARGETS})
  target_link_libraries(${ELEMENTSLIB} LINK_PRIVATE ${CPPQED_LIBRARIES} ${CPPQEDelements_LIBRARIES})

  set_target_properties(${ELEMENTSLIB} PROPERTIES
        PUBLIC_HEADER "${${PROJECT_NAME}_PUBLIC_HEADERS}"
        INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR}
        VERSION ${CPPQED_ABI_MAJOR}.${CPPQED_ABI_MINOR}.${CPPQED_ABI_MICRO}
        SOVERSION ${CPPQED_ABI_MAJOR}
  )

  install(TARGETS ${ELEMENTSLIB}
          EXPORT CPPQED${PROJECT_NAME}Targets
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${ELEMENTS_INCLUDE_SUBDIR}
          COMPONENT shlib
  )

  # Add all targets to the build-tree export set
  if(DEFINED CPPQED_MONOLITHIC)
    # workaround: use APPEND because we do not want to export ${CPPQED_LIBRARIES} here
    # otherwise external project would import this target twice
    file(REMOVE "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}Targets.cmake")
    export(TARGETS ${ELEMENTSLIB} FILE "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}Targets.cmake" APPEND)
  else(DEFINED CPPQED_MONOLITHIC)
    export(TARGETS ${ELEMENTSLIB} FILE "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}Targets.cmake")
  endif(DEFINED CPPQED_MONOLITHIC)

  option(REGISTRY "Register build trees in the cmake registry so that other projects can find them." ON)
  if(REGISTRY)
    export(PACKAGE CPPQED${PROJECT_NAME})
  endif(REGISTRY)

  # Create the CPPQEDConfig.cmake
  # ... for the build tree
  set(CONF_INCLUDE_DIRS ${PROJECT_BINARY_DIR})
  foreach(d ${ELEMENTS_SOURCE_DIRS})
    set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/${d}) 
  endforeach(d)
  set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
  configure_package_config_file(${CPPQED_CMAKE_DIR}/ElementsTemplateConfig.cmake.in "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION  "${PROJECT_BINARY_DIR}"
    PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR
  )
  write_basic_package_version_file(${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}ConfigVersion.cmake 
    VERSION ${CPPQED_VERSION_MAJOR}.${CPPQED_VERSION_MINOR}
    COMPATIBILITY ExactVersion
  )

  # ... and for the installation tree
  set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR}/${ELEMENTS_INCLUDE_SUBDIR})
  set(CONF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/${ELEMENTS_CMAKE_SUBDIR})
  configure_package_config_file(${CPPQED_CMAKE_DIR}/ElementsTemplateConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQED${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CONF_CMAKE_DIR}"
    PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR
  )

  # Install the CPPQEDConfig.cmake and CPPQEDConfigVersion.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CPPQED${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/${ELEMENTS_CMAKE_SUBDIR}" COMPONENT dev)

  # Install the export set for use with the install-tree
  install(EXPORT CPPQED${PROJECT_NAME}Targets DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/${ELEMENTS_CMAKE_SUBDIR}" COMPONENT dev)

endmacro()

#! \brief Generate C++ source code containing version information as global variables.
#! \return This function sets `CONF_GIT_SHA1` and `CONF_VERSION` (c.f. get_git_head_revision()).
#!
#! This extracts the current git commit hash by using get_git_head_revision(), and generates two files
#! `${PROJECT_NAME}_version.cc` and `${PROJECT_NAME}_version.h`. For example, with a project `elements`, clients including
#! `elements_version.h` have access to these global entities:
#!
#! - `char[] g_CPPQEDelements_GIT_SHA1`: the git commit hash value
#! - `char[] g_CPPQEDelements_VERSION[]`: the project version
#! - `string cppqed_@PROJECT_NAME@_version()`: a function returning both theses values in a human readable format.
#!
#! Generated projects (c.f. elements_project() and scripts_project()) automatically call this
#! function and link in the generated source file.
function(generate_version_files)
  get_git_head_revision(REFSPEC CONF_GIT_SHA1)
  if(CONF_GIT_SHA1 STREQUAL "GITDIR-NOTFOUND")
    set(CONF_GIT_SAH1 "built from source package")
  endif()
  set(CONF_VERSION ${CPPQED_VERSION})
  if(DEFINED CPPQED_CMAKE_DIR)
    set(_dir ${CPPQED_CMAKE_DIR})
  elseif(DEFINED CPPQED_CMAKE_MODULE_PATH)
    set(_dir ${CPPQED_CMAKE_MODULE_PATH})
  endif()
  configure_file("${_dir}/version.cc.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.cc" @ONLY)
  configure_file("${_dir}/version.h.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.h" @ONLY)
  set(CONF_GIT_SHA1 ${CONF_GIT_SHA1} PARENT_SCOPE)
  set(CONF_VERSION ${CONF_VERSION} PARENT_SCOPE)
endfunction()

#! \brief Initialize a scripts project.
#!
#! \param ELEMENTS_PROJECT `<proj_name>` (optional) Scripts depend on this custom elements project (C++QED elements project does
#!                                     not need to be supplied here). Repeat for several custom elements projects.
#!
#! Scripts which reside in the top level project directory and have a `.cc` extension are picked up for compilation. The following
#! actions are performed:
#!
#! - Find C++QED core, elements all custom element projects passed into the macro and set them up as dependencies.
#! - Generate version information for the scripts project, see generate_version_files() for details.
#! - Gather version informations from all C++QED dependencies (core, elements and custom element projects) and compile them in a
#!  generated source file `component_versions.cc`. Scripts including the header file `component_versions.h` have access to the function
#!  `std::string cppqed_component_versions()`, which returns the version information of all used components in a human readable form.
#! - Create a target for every script found. If the script (without file extension) is listed in
#!   `EXCLUDE_FROM_ALL_SCRIPTS`, this script will not be built automatically. If it is listed in `EXCLUDE_SCRIPTS`, no target will be created.
#! - Exclude all scripts listed in `NEED_FLENS` (without file extension) from compilation, if the current C++QED library does not support FLENS.
#! - Create a target `${PROJECT_NAME}_all` which compiles all scripts not listed in `EXCLUDE_FROM_ALL_SCRIPTS`.
macro(scripts_project)
  # find CPPQED elements project
  find_package(CPPQEDelements ${CPPQED_ID} REQUIRED)
  include_directories(${CPPQEDelements_INCLUDE_DIRS})
  set(ALL_ELEMENTS_LIBRARIES ${CPPQEDelements_LIBRARIES})

  include_directories(${PROJECT_BINARY_DIR})

  CPPQED_SETUP()

  generate_version_files()

  set(CONF_COMPONENT_VERSIONS "cppqed_core_version()+cppqed_elements_version()")
  set(CONF_COMPONENT_VERSIONS_INCLUDES "#include \"core_version.h\"\n#include \"elements_version.h\"\n")

  # find additional elements projects
  foreach(elements ${ARGV})
    find_package(CPPQED${elements} ${CPPQED_ID} REQUIRED)
    set(ALL_ELEMENTS_LIBRARIES ${ALL_ELEMENTS_LIBRARIES} ${CPPQED${elements}_LIBRARIES})
    include_directories(${CPPQED${elements}_INCLUDE_DIRS})
    set(CONF_COMPONENT_VERSIONS_INCLUDES "${CONF_COMPONENT_VERSIONS_INCLUDES}#include \"${elements}_version.h\"\n")
    set(CONF_COMPONENT_VERSIONS "${CONF_COMPONENT_VERSIONS}+cppqed_${elements}_version()")
  endforeach()

  set(CONF_COMPONENT_VERSIONS_INCLUDES "${CONF_COMPONENT_VERSIONS_INCLUDES}#include \"${PROJECT_NAME}_version.h\"")
  set(CONF_COMPONENT_VERSIONS "${CONF_COMPONENT_VERSIONS}+cppqed_${PROJECT_NAME}_version()")

  configure_file(${CPPQED_CMAKE_DIR}/component_versions.cc.in ${PROJECT_BINARY_DIR}/component_versions.cc @ONLY)
  configure_file(${CPPQED_CMAKE_DIR}/component_versions.h.in ${PROJECT_BINARY_DIR}/component_versions.h @ONLY)

  add_library(${PROJECT_NAME}_versions_obj OBJECT ${PROJECT_BINARY_DIR}/component_versions.cc ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_version.cc)
  set_target_properties(${PROJECT_NAME}_versions_obj PROPERTIES POSITION_INDEPENDENT_CODE On)

  # create all scripts targets
  file(GLOB SCRIPTS . *.cc)
  set(EXCLUDE_FROM_ALL_SCRIPTS ${EXCLUDE_FROM_ALL_SCRIPTS} ${CPPQED_EXCLUDE_SCRIPTS})
  foreach(s ${SCRIPTS})
    get_filename_component(SCRIPT ${s} NAME_WE)
    list(FIND NEED_FLENS ${SCRIPT} F)
    list(FIND EXCLUDE_SCRIPTS ${SCRIPT} EX)
    list(FIND EXCLUDE_FROM_ALL_SCRIPTS ${SCRIPT} NOT_ALL)
    set(exclude Off)
    if(NOT EX EQUAL -1)
      set(exclude On)
    endif()
    if(NOT (CPPQED_FLENS_FOUND OR F EQUAL -1))
      set(exclude On)
    endif()
    if(NOT exclude)
      add_executable(${SCRIPT} ${s} $<TARGET_OBJECTS:${PROJECT_NAME}_versions_obj>)
      target_link_libraries(${SCRIPT} ${CPPQED_LIBRARIES} ${ALL_ELEMENTS_LIBRARIES})
      set_target_properties(${SCRIPT} PROPERTIES DEBUG_POSTFIX _d)
      if(NOT_ALL EQUAL -1)
        set(SCRIPTNAMES ${SCRIPTNAMES} ${SCRIPT})
        install(TARGETS ${SCRIPT}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
      else()
        set_target_properties(${SCRIPT} PROPERTIES EXCLUDE_FROM_ALL ON)
      endif()
    endif()
  endforeach(s)
  # add target "scripts"
  add_custom_target(${PROJECT_NAME}_all)
  if(SCRIPTNAMES)
    add_dependencies(${PROJECT_NAME}_all ${SCRIPTNAMES})
  endif()
endmacro()

#! \brief Generate a doxygen documentation target for the current project.
#!
#! \param target_prefix The name of the documentation target will be ${target_prefix}_doc
#! \param tagfiles A list with full paths to doxygen tagfiles of other projects.
#! \param dependencies (optional) All remaining arguments will be treated as targets
#!  this document target should depend on.
#!
#! If doxygen or dot is not found on the system, this macro does nothing.
#! Otherwise it sets the variable `TAGFILES` to a list which can be used in the TAGFILES option of a
#! Doxyfile. The HTML location corresponding to each tagfile is expected in a html subdirectory of the
#! directory where the tagfile resides. All paths will be converted to relative paths so that the
#! resulting documentation can be relocated.
#!
#! The template `doc/Doxyfile` will be copied to `${PROJECT_NAME}_DOC_DIR`, expanding all @-variables within.
macro(cppqed_documentation target_prefix tagfiles)

  set(tagfiles ${tagfiles})

  if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    set(doc_depends ${ARGN})
    set(DOXYGEN_HEADER_FILE ${cppqed_SOURCE_DIR}/doc/header.html)
    set(DOXYGEN_CSS_FILE ${cppqed_SOURCE_DIR}/doc/stylesheet.css)
    file(MAKE_DIRECTORY ${${PROJECT_NAME}_DOC_DIR})
    while(tagfiles)
      list(GET tagfiles 0 tagfile)
      list(REMOVE_AT tagfiles 0)
      file(RELATIVE_PATH relative_tagfile ${${PROJECT_NAME}_DOC_DIR} ${tagfile})
      get_filename_component(relative_location ${relative_tagfile} PATH)
      set(TAGFILES "${TAGFILES} ${relative_tagfile}=../${relative_location}/html")
    endwhile()
    set(local_doc_dir ${${PROJECT_NAME}_DOC_DIR})
    configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${${PROJECT_NAME}_DOC_DIR}/Doxyfile @ONLY)
    add_custom_target(${target_prefix}doc ${DOXYGEN_EXECUTABLE} ${${PROJECT_NAME}_DOC_DIR}/Doxyfile
      WORKING_DIRECTORY ${${PROJECT_NAME}_DOC_DIR}
      COMMENT "Generating APIs documentation with Doxygen"
      DEPENDS ${doc_depends}
    )
    install(DIRECTORY ${${PROJECT_NAME}_DOC_DIR}/html
          DESTINATION ${CPPQED_DOC_DIR}/${PROJECT_NAME}
          OPTIONAL
    )
  endif()
endmacro()

#! @}
