testdir(physics fewer_scripts)

# Scripts needed for PTLA testcase

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
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PTLA_Evolved;ph_PTLA_CPPQED;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch;ph_PLQ_CPPQEDuip;ph_PLQ_CPPQEDsch;ph_PLQ_CPPQED")

if(NOT XCODE)
  add_custom_target(check_PTLA COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PTLA_Evolved\\|ph_PTLA_CPPQED\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch\\|ph_PLQ_CPPQEDuip\\|ph_PLQ_CPPQEDsch\\|ph_PLQ_CPPQED\\|compare_PTLA" DEPENDS ${TESTSDEPEND})
endif()

# Scripts needed for PTLAHL testcase

declaretest(ph_PTLAHL_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_EvolvedHL>)

declaretest(ph_PTLAHL_Si)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(compare_PTLAHL)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PTLAHL_Evolved;ph_PTLAHL_Si")
if(NOT XCODE)
  add_custom_target(check_PTLAHL COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PTLAHL_Evolved\\|ph_PTLAHL_Si\\|compare_PTLAHL")
endif()

# Scripts needed for PLM testcase

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
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_CPPQED;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")
if(NOT XCODE)
  add_custom_target(check_PLM COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_CPPQED\\|ph_PLM_CPPQEDsch\\|ph_PLM_CPPQEDuip\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch\\|compare_PLM")
endif()

declaretest(compare_PLM_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_Evolved;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")
if(NOT XCODE)
  add_custom_target(check_PLM_Ev COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_Evolved\\|ph_PLM_CPPQEDsch\\|ph_PLM_CPPQEDuip\\|ph_QbitMode_CPPQED\\|ph_QbitMode_CPPQEDsch")
endif()

# Scripts needed for Decay testcase

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
if(NOT XCODE)
  add_custom_target(check_DecayQB COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_Decay\\|ph_DecayFock\\|compare_DecayQB")
endif()

declaretest(compare_DecayMode)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Decay;ph_DecayFock")
if(NOT XCODE)
  add_custom_target(check_DecayMode COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_Decay\\|ph_DecayFock\\|compare_DecayMode")
endif()

# Scripts needed for QMJ_Int testcase

declaretest(ph_QMJ_Int_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_Evolved>)

declaretest(ph_QMJ_Int_Matrix)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_Matrix>)

declaretest(ph_QMJ_Int)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QMJ_IntSch)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QMJ_IntUIP)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(compare_QMJ_Int)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_QMJ_Int_Matrix;ph_QMJ_Int;ph_QMJ_IntSch;ph_QMJ_IntUIP")
if(NOT XCODE)
  add_custom_target(check_QMJ_Int COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_QMJ_Int_Matrix\\|ph_QMJ_Int\\|ph_QMJ_IntSch\\|ph_QMJ_IntUIP\\|compare_QMJ_Int\\|compare_QMJ_Int")
endif()

declaretest(compare_QMJ_Int_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_QMJ_Int_Ev;ph_QMJ_Int;ph_QMJ_IntSch")
if(NOT XCODE)
  add_custom_target(check_QMJ_Int_Ev COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_QMJ_Int_Ev\\|ph_QMJ_Int\\|ph_QMJ_IntSch\\|compare_QMJ_Int_Ev")
endif()

# Scripts needed for PLM_FT testcase

declaretest(ph_PLM_FT)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:PumpedLossyMode_C++QED>)

declaretest(compare_PLM_FT)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_FT")
if(NOT XCODE)
  add_custom_target(check_PLM_FT COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_FT\\|${TESTNAME}")
endif()

declaretest(compare_PLM_Rho)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=StateComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_FT")
if(NOT XCODE)
  add_custom_target(check_PLM_Rho COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_PLM_FT\\|${TESTNAME}")
endif()

# Scripts needed for FreePart testcase

declaretest(ph_FreeParticle)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:SingleParticle>)

declaretest(ph_FreeParticleSch)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:SingleParticle>)

declaretest(compare_Free_X)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_FreeParticle;ph_FreeParticleSch")

declaretest(compare_Free_VarX)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_FreeParticle;ph_FreeParticleSch")
if(NOT XCODE)
  add_custom_target(check_Free COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_FreeParticle\\|ph_FreeParticleSch\\|compare_Free_X\\|compare_Free_VarX")
endif()

# Scripts needed for Ring testcase

declaretest(ph_Ring)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_Ring_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(ph_RingSC)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_RingSC_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(ph_RingS2S)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_RingS2S_Ev)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(compare_Ring)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Ring_Ev;ph_Ring")
declaretest(compare_RingSC)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_RingSC_Ev;ph_RingSC")
declaretest(compare_RingS2S)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_RingS2S_Ev;ph_RingS2S")

if(NOT XCODE)
  add_custom_target(check_Ring COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  -R "ph_Ring\\|ph_Ring_Ev\\|ph_RingSC\\|ph_RingSC_EV\\|ph_RingS2S\\|ph_RingS2S_Ev\\compare_Ring\\|compare_RingSC\\|compare_RingS2S")
endif()

