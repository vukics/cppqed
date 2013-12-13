project(testing)

include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py)

add_test(NAME testdriver
         WORKING_DIRECTORY ${testing_BINARY_DIR}
         COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT} )