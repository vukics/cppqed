file(GLOB SCRIPTS *.cc)

if ( NOT CPPQED_HAS_FLENS )
  # this cannot be built without FLENS support
  list(REMOVE_ITEM SCRIPTS ${CMAKE_SOURCE_DIR}/src/GeneralDickeImaginaryEvolution.cc)
endif( NOT CPPQED_HAS_FLENS )

foreach(s ${SCRIPTS})
  get_filename_component(SCRIPT ${s} NAME_WE)
  set(SCRIPTNAMES ${SCRIPTNAMES} ${SCRIPT})
  add_executable(${SCRIPT} ${s})
  target_link_libraries(${SCRIPT} ${CPPQED_LIBRARIES})
  #install(TARGETS ${SCRIPT}
  #        RUNTIME DESTINATION bin)
endforeach(s)

# add target "scripts"
add_custom_target(scripts)
add_dependencies(scripts ${SCRIPTNAMES})

# add target "fewer_scripts"
add_custom_target(fewer_scripts)
add_dependencies(fewer_scripts PTLA_Evolved PTLA_C++QED PumpedLossyQbit PTLA_EvolvedHL PumpedLossyMode_Evolved PumpedLossyMode_C++QED QbitMode_C++QED QbitMode_Evolved QbitMode_Matrix)


