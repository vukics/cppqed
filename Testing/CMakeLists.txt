project(testing)

include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

add_custom_target(check
                  DEPENDS check_run testdriver.conf.in
)

set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf)

# Testing the testdriver itself

add_test(NAME testdriver
         COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT} )

add_subdirectory(run)