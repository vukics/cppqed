aux_source_directory(. SCRIPTS)

include_directories(${ELEMENTS_INCLUDE_DIRS})

if ( NOT flens_FOUND )
  # this cannot be built without FLENS support
  list(REMOVE_ITEM SCRIPTS ./GeneralDickeImaginaryEvolution.cc)
endif( NOT flens_FOUND )

# checks if libutils should be linked with FLENS
if ( CMAKE_CXX_COMPILER_VERSION MATCHES " 4.2" )
  message(STATUS "Detected gcc version 4.2, removing some scripts from build targets.")
  list(REMOVE_ITEM SCRIPTS ./2particles1mode.cc ./Ca40InCavity.cc ./QbitMode_Matrix.cc ./MLJC.cc )
endif ( CMAKE_CXX_COMPILER_VERSION MATCHES " 4.2" )

foreach(s ${SCRIPTS})
  get_filename_component(SCRIPT ${s} NAME_WE)
  add_executable(${SCRIPT} ${s})
  target_link_libraries(${SCRIPT} el)
  install(TARGETS ${SCRIPT}
	  RUNTIME DESTINATION bin)
endforeach(s)

