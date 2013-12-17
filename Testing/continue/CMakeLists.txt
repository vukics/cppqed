add_custom_target(check_continue
                  COMMAND ${CMAKE_CTEST_COMMAND} -V
                  DEPENDS 1particle1mode
)

set(TESTNAME continue_1particle1mode_ndt)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Continuer --script=$<TARGET_FILE:1particle1mode>)