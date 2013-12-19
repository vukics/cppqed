project(testing)

# Please avoid any special characters in test names (+,^,.,etc.)
# The ctest command might interpret them as regular expressions


include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

add_custom_target(check DEPENDS io check_run check_continue check_boost testdriver.conf.in)

macro(testdir dir)
  set(TESTSDEPEND ${ARGN})
  add_custom_target(check_${dir}
                  COMMAND ${CMAKE_CTEST_COMMAND} -V
                  DEPENDS ${TESTSDEPEND}
  )
  set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf)
endmacro()

macro(declaretest name)
  set(TESTNAME ${name})
  add_custom_target(${name} COMMAND ${CMAKE_CTEST_COMMAND} -V -R "^${name}$" DEPENDS ${TESTSDEPEND})
endmacro()

add_subdirectory(run)
add_subdirectory(continue)
add_subdirectory(compilefail)
add_subdirectory(boost)