testdir(cpypyqed)

set(ONDEMAND ${PROJECT_BINARY_DIR}/ondemand_temp)

configure_file(config.txt.in ${ONDEMAND}/config.txt @ONLY)

declaretest(cpypyqed_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/2particles1mode.py
            --cpypyqed_builddir=${ONDEMAND})

