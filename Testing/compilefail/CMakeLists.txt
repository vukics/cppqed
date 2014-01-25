testdir(compilefail)

CPPQED_SETUP()

add_library(CompositeError2 OBJECT EXCLUDE_FROM_ALL CompositeError2.cc)
declaretest(fail_CompositeError2)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=CompositeError2)
