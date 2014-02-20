testdir(physics PTLA_Evolved PTLA_C++QED QbitMode_C++QED PumpedLossyQbit PTLA_EvolvedHL
                PumpedLossyMode_C++QED PumpedLossyMode_Evolved)

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

# Scripts needed for PLM comparison

declaretest(ph_PLM_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyMode_Evolved>)

declaretest(ph_PLM_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyMode_C++QED>)

declaretest(ph_PLM_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyMode_C++QED>)

declaretest(ph_PLM_CPPQEDuip)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyMode_C++QED>)

declaretest(compare_PLM)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Comparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_CPPQED;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")
add_custom_target(check_PLM COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_CPPQED\\|ph_PLM_CPPQEDsch\\|ph_PLM_CPPQEDuip\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch\\|compare_PLM")

declaretest(compare_PLM_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Comparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_Evolved;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")
add_custom_target(check_PLM_Ev COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_Evolved\\|ph_PLM_CPPQEDsch\\|ph_PLM_CPPQEDuip\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch")


# Scripts needed for Decay comparison

declaretest(ph_Decay)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_DecayFock)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(compare_DecayQB)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Decay;ph_DecayFock")
add_custom_target(check_DecayQB COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_Decay\\|ph_DecayFock\\|compare_DecayQB")

declaretest(compare_DecayMode)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Decay;ph_DecayFock")
add_custom_target(check_DecayMode COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_Decay\\|ph_DecayFock\\|compare_DecayMode")
