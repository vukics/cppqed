# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

declaretest(ph_DTLA_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DTLA_Evolved>)
set_tests_properties(${TESTNAME} PROPERTIES LABELS "Physics")

declaretest(ph_DTLA_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DTLA_C++QED>)

declaretest(ph_QbitMode_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QbitMode_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_PLQ_CPPQEDuip)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeQbit>)

declaretest(ph_PLQ_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeQbit>)

declaretest(ph_PLQ_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeQbit>)

declaretest(compare_DTLA)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_DTLA_Evolved;ph_DTLA_CPPQED;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch;ph_PLQ_CPPQEDuip;ph_PLQ_CPPQEDsch;ph_PLQ_CPPQED")

# Scripts needed for DTLAHL testcase

declaretest(ph_DTLAHL_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DTLA_EvolvedHL>)

declaretest(ph_DTLAHL_Si)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DTLA_C++QED>)

declaretest(compare_DTLAHL)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "ph_DTLAHL_Evolved;ph_DTLAHL_Si")

# Scripts needed for PLM testcase

declaretest(ph_PLM_Evolved)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeMode_Evolved>)

declaretest(ph_PLM_CPPQED)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeMode_C++QED>)

declaretest(ph_PLM_CPPQEDsch)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeMode_C++QED>)

declaretest(ph_PLM_CPPQEDuip)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeMode_C++QED>)

declaretest(compare_PLM)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_CPPQED;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")


declaretest(compare_PLM_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_Evolved;ph_PLM_CPPQEDsch;ph_PLM_CPPQEDuip;ph_QbitMode_CPPQED;ph_QbitMode_CPPQEDsch")


# Scripts needed for Decay testcase

declaretest(ph_Decay)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_DecayFock)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(compare_DecayQB)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Decay;ph_DecayFock")

declaretest(compare_DecayMode)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Decay;ph_DecayFock")

# Scripts needed for QMJ_Int testcase

declaretest(ph_QMJ_Int_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_Evolved>)

declaretest(ph_QMJ_Int_Matrix)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_Matrix>)

declaretest(ph_QMJ_Int)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QMJ_IntSch)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(ph_QMJ_IntUIP)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:QbitMode_C++QED>)

declaretest(compare_QMJ_Int)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_QMJ_Int_Matrix;ph_QMJ_Int;ph_QMJ_IntSch;ph_QMJ_IntUIP")

declaretest(compare_QMJ_Int_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_QMJ_Int_Ev;ph_QMJ_Int;ph_QMJ_IntSch")

# Scripts needed for PLM_FT testcase

declaretest(ph_PLM_FT)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:DrivenDissipativeMode_C++QED>)

declaretest(compare_PLM_FT)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_FT")

declaretest(compare_PLM_Rho)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=StateComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_PLM_FT")

# Scripts needed for FreePart testcase

declaretest(ph_FreeParticle)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:SingleParticle>)

declaretest(ph_FreeParticleSch)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:SingleParticle>)

declaretest(compare_Free_X)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_FreeParticle;ph_FreeParticleSch")

declaretest(compare_Free_VarX)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=FunctionComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_FreeParticle;ph_FreeParticleSch")

# Scripts needed for Ring testcase

declaretest(ph_Ring)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_Ring_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(ph_RingSC)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_RingSC_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(ph_RingS2S)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring>)

declaretest(ph_RingS2S_Ev)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:Ring_Evolved>)

declaretest(compare_Ring)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_Ring_Ev;ph_Ring")
  
declaretest(compare_RingSC)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_RingSC_Ev;ph_RingSC")

declaretest(compare_RingS2S)
add_test( NAME ${TESTNAME} COMMAND ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=TrajectoryComparer)
set_tests_properties(${TESTNAME} PROPERTIES
  DEPENDS "ph_RingS2S_Ev;ph_RingS2S")

