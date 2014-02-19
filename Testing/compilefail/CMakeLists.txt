testdir(compilefail)

CPPQED_SETUP()

file(GLOB SOURCEFILES . *.cc)

foreach(SOURCEFILE ${SOURCEFILES})
  get_filename_component(S ${SOURCEFILE} NAME_WE)
  add_library( ${S} OBJECT EXCLUDE_FROM_ALL ${SOURCEFILE} )
  declaretest( fail_${S} )
  add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
              --test=${TESTNAME} --testclass=CompileTarget --script=${S} )
endforeach()
