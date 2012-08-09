aux_source_directory(. SCRIPTS)

include_directories(${ELEMENTS_INCLUDE_DIRS})

IF( NOT flens_FOUND )
  # this cannot be built without FLENS support
  LIST(REMOVE_ITEM SCRIPTS ./GeneralDickeImaginaryEvolution.cc)
ENDIF( NOT flens_FOUND )

foreach(s ${SCRIPTS})
  get_filename_component(SCRIPT ${s} NAME_WE)
  add_executable(${SCRIPT} ${s})
  target_link_libraries(${SCRIPT} el)
endforeach(s)