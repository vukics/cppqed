testdir(run 1particle1mode PTLA_C++QED PTLA_Evolved 3qbits)

# PTLA_C++QED as a representative of simulating a free

declaretest(run_PTLA_CPPQED_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(run_PTLA_CPPQED_dt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(run_PTLA_CPPQED_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_C++QED>)

declaretest(run_PTLA_Evolved_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_Evolved>)

declaretest(run_PTLA_Evolved_dt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_Evolved>)

declaretest(run_PTLA_Evolved_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:PTLA_Evolved>)

# 1particle1mode as a representative of simulating a binary interaction

declaretest(run_1particle1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)

declaretest(run_1particle1mode_dt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)

declaretest(run_1particle1mode_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)

# 3qbits as a representative of simulating a composite

declaretest(run_3qbits_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:3qbits>)

declaretest(run_3qbits_dt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:3qbits>)

declaretest(run_3qbits_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:3qbits>)
