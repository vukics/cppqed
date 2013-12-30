testdir(cpypyqed 2particles1mode)

set(ONDEMAND ${PROJECT_BINARY_DIR}/ondemand_temp)

configure_file(config.txt.in ${ONDEMAND}/config.txt @ONLY)

declaretest(cpypyqed_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonContinuer --script=${cpypyqed_SOURCE_DIR}/scripts/2particles1mode.py
            --cpypyqed_builddir=${ONDEMAND})

declaretest(cppqed_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Continuer --script=$<TARGET_FILE:2particles1mode>)

declaretest(verify_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Verifier)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cppqed_2particles1mode_dc;cpypyqed_2particles1mode_dc")