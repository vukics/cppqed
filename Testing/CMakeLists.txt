project(testing)

# Please avoid any special characters in test names (+,^,.,etc.)
# The ctest command might interpret them as regular expressions


include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2.7 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

add_custom_target(check_fast DEPENDS check_compilefail check_run check_continue check_boost testdriver.conf.in)
add_custom_target(check DEPENDS check_fast check_cpypyqed)

set(TEST_PARALLEL 1 CACHE STRING "Number of parallel test jobs for the test suite.")
set_property(CACHE TEST_PARALLEL PROPERTY STRINGS 1 2 3 4 5 6 7 8 9 10)
set(CTEST_J "-j${TEST_PARALLEL}")
option(TEST_VERBOSE "Verbose output in test suite." ON)
if(TEST_VERBOSE)
  set(CTEST_V "-V")
endif()

macro(testdir dir)
  set(TESTSDEPEND cpypyqed ${ARGN})
  add_custom_target(check_${dir}
                  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  DEPENDS ${TESTSDEPEND}
  )
  set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --configuration=$<CONFIGURATION>)
endmacro()

macro(declaretest name)
  set(TESTNAME ${name})
  add_custom_target(${name} COMMAND ${CMAKE_CTEST_COMMAND} -V -R "^${name}$" DEPENDS ${TESTSDEPEND})
endmacro()

add_subdirectory(run)
add_subdirectory(continue)
add_subdirectory(compilefail)
add_subdirectory(boost)
add_subdirectory(cpypyqed)
