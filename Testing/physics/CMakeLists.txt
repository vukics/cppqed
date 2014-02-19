testdir(physics PTLA_Evolved PTLA_C++QED QbitMode_C++QED)

declaretest(ph_PTLA_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_Evolved>)

declaretest(ph_PTLA_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(ph_QbitMode_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QbitMode_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(compare_PTLA)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Comparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PTLA_Evolved;ph_PTLA_CPPQED;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")
