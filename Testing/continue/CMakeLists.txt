# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

declaretest(continue_1particle1mode_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Continuer --script=$<TARGET_FILE:1particle1mode>)

declaretest(verify_continue_1particle1mode_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS continue_1particle1mode_ndt)
