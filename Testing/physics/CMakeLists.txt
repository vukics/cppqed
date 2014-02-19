testdir(physics PTLA_Evolved PTLA_C++QED QbitMode_C++QED PumpedLossyQbit PTLA_EvolvedHL)

# Scripts needed for PTLA comparison

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

declaretest(ph_PLQ_CPPQEDuip)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyQbit>)

declaretest(ph_PLQ_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyQbit>)

declaretest(ph_PLQ_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyQbit>)

declaretest(compare_PTLA)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Comparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PTLA_Evolved;ph_PTLA_CPPQED;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch;ph_PLQ_CPPQEDuip;ph_PLQ_CPPQEDsch;ph_PLQ_CPPQED")

add_custom_target(check_PTLA COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PTLA_Evolved\\|ph_PTLA_CPPQED\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch\\|ph_PLQ_CPPQEDuip\\|ph_PLQ_CPPQEDsch\\|ph_PLQ_CPPQED\\|compare_PTLA" DEPENDS ${TESTSDEPEND})

# Scripts needed for PTLAHL comparison

declaretest(ph_PTLAHL_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_EvolvedHL>)

declaretest(ph_PTLAHL_Si)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(compare_PTLAHL)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Comparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PTLAHL_Evolved;ph_PTLAHL_Si")
add_custom_target(check_PTLAHL COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PTLAHL_Evolved\\|ph_PTLAHL_Si\\|compare_PTLAHL")
