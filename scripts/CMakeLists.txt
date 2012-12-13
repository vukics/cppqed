file(GLOB SCRIPTS *.cc)
include_directories(${ELEMENTS_INCLUDE_DIRS})

if ( NOT ENABLE_FLENS )
  # this cannot be built without FLENS support
  list(REMOVE_ITEM SCRIPTS ${CMAKE_SOURCE_DIR}/scripts/GeneralDickeImaginaryEvolution.cc)
endif( NOT ENABLE_FLENS )

# remove 6qbits which causes problems on some compilers
list(REMOVE_ITEM SCRIPTS ${CMAKE_SOURCE_DIR}/scripts/6qbits.cc)

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
  add_executable(${SCRIPT} ${s})
  target_link_libraries(${SCRIPT} ${CPPQEDLIB})
  install(TARGETS ${SCRIPT}
	  RUNTIME DESTINATION bin)
endforeach(s)

