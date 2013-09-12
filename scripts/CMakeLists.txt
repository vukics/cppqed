file(GLOB SCRIPTS *.cc)

include_directories(${ELEMENTS_INCLUDE_DIRS})

if ( NOT ENABLE_FLENS )
  # this cannot be built without FLENS support
  list(REMOVE_ITEM SCRIPTS ${CMAKE_SOURCE_DIR}/scripts/GeneralDickeImaginaryEvolution.cc)
endif( NOT ENABLE_FLENS )

# checks if libutils should be linked with FLENS
if ( CMAKE_CXX_COMPILER_VERSION MATCHES " 4.2" )
  message(STATUS "Detected gcc version 4.2, removing some scripts from build targets.")
  list(REMOVE_ITEM SCRIPTS ${CMAKE_SOURCE_DIR}/scripts/2particles1mode.cc
      ${CMAKE_SOURCE_DIR}/scripts/Ca40InCavity.cc
      ${CMAKE_SOURCE_DIR}/scripts/QbitMode_Matrix.cc
      ${CMAKE_SOURCE_DIR}/scripts/MLJC.cc )
endif ( CMAKE_CXX_COMPILER_VERSION MATCHES " 4.2" )

foreach(s ${SCRIPTS})
  get_filename_component(SCRIPT ${s} NAME_WE)
  set(SCRIPTNAMES ${SCRIPTNAMES} ${SCRIPT})
  add_executable(${SCRIPT}  EXCLUDE_FROM_ALL ${s})
  target_link_libraries(${SCRIPT} ${CPPQEDLIB})
  #install(TARGETS ${SCRIPT}
  #        RUNTIME DESTINATION bin)
endforeach(s)

# add target "scripts"
add_custom_target(scripts)
add_dependencies(scripts ${SCRIPTNAMES})

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix)


