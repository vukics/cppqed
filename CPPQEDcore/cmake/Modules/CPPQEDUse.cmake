
MACRO(CPPQED_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-ignored-qualifiers -Wno-sign-compare -Wno-overloaded-virtual")
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-local-type-template-args")
  endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
ENDMACRO(CPPQED_CXX_FLAGS)

# This function populates the variable ${namespace}_PUBLIC_HEADERS with all the
# headers in the list ${source_dirs}
function(gather_includes namespace source_dirs)
  foreach(d ${${source_dirs}})
    file(GLOB INC ${d}/*.h ${d}/*.tcc)
    set(${namespace}_PUBLIC_HEADERS ${${namespace}_PUBLIC_HEADERS} ${INC})
  endforeach(d)
  set(${namespace}_PUBLIC_HEADERS ${${namespace}_PUBLIC_HEADERS} PARENT_SCOPE)
endfunction()

# This macro creates a target ${PROJECT_NAME}_${name}_objs where ${name} is the current
# subdirectory and populates it with all source files in this directory. It then adds this 
# target to the parent scope variable ${OBJ_TARGETS}.
macro(create_object_target)
  get_filename_component(name ${CMAKE_CURRENT_LIST_DIR} NAME)
  include_directories(.)
  foreach(d in ${${name}_NEEDS})
    include_directories(../${d})
  endforeach(d)
  aux_source_directory(. ${name}_srcs)
  add_library(${PROJECT_NAME}_${name}_objs OBJECT ${${name}_srcs})
  set_target_properties(${PROJECT_NAME}_${name}_objs PROPERTIES POSITION_INDEPENDENT_CODE On)
  set(OBJ_TARGETS ${OBJ_TARGETS} "\$<TARGET_OBJECTS:${PROJECT_NAME}_${name}_objs>" PARENT_SCOPE)
endmacro()

macro(CPPQED_SETUP)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  set(CMAKE_DEBUG_POSTFIX "_d")

  cppqed_cxx_flags()

  set(${CMAKE_CXX_FLAGS_DEBUG} "${CMAKE_CXX_FLAGS_DEBUG} -DBZ_DEBUG")
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

macro(elements_project)
  if(${PROJECT_NAME} STREQUAL elements)
    if(NOT DEFINED ORIGINAL_ELEMENTS_PROJECT)
      message(FATAL_ERROR "The project cannot be named 'elements', as this is the name of the CPPQED elements project.")
    endif()
  else()
    find_package(CPPQEDelements ${CPPQED_ID} REQUIRED)
    include_directories(${CPPQEDelements_INCLUDE_DIRS})
  endif()

  CPPQED_SETUP()

  foreach(d ${ELEMENTS_SOURCE_DIRS})
    add_subdirectory(${d})
  endforeach(d)

  set(ELEMENTS_CMAKE_SUBDIR "cmake/CPPQED${PROJECT_NAME}-${CPPQED_ID}")
  set(ELEMENTS_INCLUDE_SUBDIR "CPPQED-${CPPQED_ID}/${PROJECT_NAME}")

  gather_includes(ELEMENTS ELEMENTS_SOURCE_DIRS)

  set(ELEMENTSLIB C++QED${PROJECT_NAME}-${CPPQED_ID})
  add_library(${ELEMENTSLIB} SHARED ${ELEMENTS_PUBLIC_HEADERS} ${OBJ_TARGETS})
  target_link_libraries(${ELEMENTSLIB} LINK_PRIVATE ${CPPQED_LIBRARIES} ${CPPQEDELEMENTS_LIBRARIES})

  set_target_properties(${ELEMENTSLIB} PROPERTIES
        PUBLIC_HEADER "${ELEMENTS_PUBLIC_HEADERS}"
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
  export(PACKAGE CPPQED${PROJECT_NAME})

  # Create the CPPQEDConfig.cmake
  # ... for the build tree
  foreach(d ${ELEMENTS_SOURCE_DIRS})
    set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/${d}) 
  endforeach(d)
  set(CONF_CMAKE_DIR ${PROJECT_BINARY_DIR})
  configure_package_config_file(${CPPQED_CMAKE_DIR}/ElementsTemplateConfig.cmake.in "${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION  "${PROJECT_BINARY_DIR}"
    PATH_VARS CONF_INCLUDE_DIRS CPPQED_THIRDPARTY_INCLUDE_DIRS CONF_CMAKE_DIR
  )
  write_basic_package_version_file(${PROJECT_BINARY_DIR}/CPPQED${PROJECT_NAME}ConfigVersion.cmake 
    VERSION ${CPPQED_MAJOR_VERSION}.${CPPQED_MINOR_VERSION}
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

macro(scripts_project)
  # find CPPQED elements project
  if(NOT DEFINED CPPQED_MONOLITHIC)
    find_package(CPPQEDelements ${CPPQED_ID} REQUIRED)
  else(NOT DEFINED CPPQED_MONOLITHIC)
    include(${CPPQED_BINARY_DIR}/CPPQEDelements/CPPQEDelementsConfig.cmake)
  endif(NOT DEFINED CPPQED_MONOLITHIC)
  include_directories(${CPPQEDelements_INCLUDE_DIRS})
  set(ALL_ELEMENTS_LIBRARIES ${CPPQEDelements_LIBRARIES})

  CPPQED_SETUP()

  # find additional elements projects
  foreach(elements ${ARGV})
    find_package(CPPQED${elements} ${CPPQED_ID} REQUIRED)
    set(ALL_ELEMENTS_LIBRARIES ${ALL_ELEMENTS_LIBRARIES} ${CPPQED${elements}_LIBRARIES})
    include_directories(${CPPQED${elements}_INCLUDE_DIRS})
  endforeach()

  # create all scripts targets
  file(GLOB SCRIPTS . *.cc)
  foreach(s ${SCRIPTS})
    get_filename_component(SCRIPT ${s} NAME_WE)
    list(FIND NEED_FLENS ${SCRIPT} F)
    list(FIND EXCLUDE_SCRIPTS ${SCRIPT} EX)
    if(EX EQUAL -1 AND (CPPQED_HAS_FLENS OR F EQUAL -1))
      set(SCRIPTNAMES ${SCRIPTNAMES} ${SCRIPT})
      add_executable(${SCRIPT} ${s})
      target_link_libraries(${SCRIPT} ${CPPQED_LIBRARIES} ${ALL_ELEMENTS_LIBRARIES})
      install(TARGETS ${SCRIPT}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
  endforeach(s)
  # add target "scripts"
  add_custom_target(scripts)
  add_dependencies(scripts ${SCRIPTNAMES})
endmacro()