project(testing)

include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py)

add_test(NAME testdriver
         COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT} )

set(TESTNAME run_1particle1mode_dc)
add_test( NAME ${TESTNAME}
          COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Runner --script=$<TARGET_FILE:1particle1mode>
)