testdir(compile)

declaretest(test_CustomElementsExample)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=C++QEDelements_custom_example-2.99)
declaretest(test_CustomScriptsExample)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=scripts_custom_example_all)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS test_CustomElementsExample)
declaretest(test_scripts_all_target)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=scripts_all)