#! \ingroup Main
#! \file
#! \brief %CMake file for the C++QED testing suite.
#!
#! The file has the following structure:

project(testing)

# Please avoid any special characters in test names (+,^,.,etc.)
# The ctest command might interpret them as regular expressions


#! \file
#! <!--#########################################################-->
#! ### Configuration of the test suite
#! <!--#########################################################-->
#!
#! Find the Python interpreter and add the targets `check_fast` and `ckeck` (c.f. \ref testsuite_running).
#! Also defines helper macros.

include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
include(${CPPQED_USE})

find_package(PythonInterp 2.7 REQUIRED)

configure_file(testdriver.conf.in testdriver.conf)

add_custom_target(check_fast DEPENDS check_compilefail check_run check_continue check_boost testdriver.conf.in)
add_custom_target(check DEPENDS check_fast check_compile check_cpypyqed)

#! \name Project options
#! @{

#! Number of parallel test jobs for the test suite, default 1.
set(TEST_PARALLEL 1 CACHE STRING "Number of parallel test jobs for the test suite.")
set_property(CACHE TEST_PARALLEL PROPERTY STRINGS 1 2 3 4 5 6 7 8 9 10)
set(CTEST_J "-j${TEST_PARALLEL}")
#! Switch verbose test output.
option(TEST_VERBOSE "Verbose output in test suite." ON)
#! @}
if(TEST_VERBOSE)
  set(CTEST_V "-V")
endif()

#! \brief Initialize a sub-directory of the test suite.
#! \ingroup Helpers
#! \param dir The sub-directory name.
#! \return `TESTSCRIPT` - set to the command by which to call the Python test driver.
#!
#! This macro adds the target `check_${dir}` and makes it depend on any target passed as argument after `dir`.
#! This target runs all tests in the current sub-directory.
#! It also sets the variable `TESTSCRIPT`, which holds the full path to the python test driver and already
#! adds some command line parameters: the configuration files and the build configuration.
macro(testdir dir)
  set(TESTSDEPEND cpypyqed ${ARGN})
  add_custom_target(check_${dir}
                  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                  DEPENDS ${TESTSDEPEND}
  )
  set(TESTSCRIPT ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --configuration=$<CONFIGURATION>)
endmacro()

#! \brief Declare a new test.
#! \ingroup Helpers
#! \param name Name of the test. Avoid any special characters (+,\,$ etc), which might be mis-interpreted in a regular expression.
#! \return `TESTNAME` - holds the name of the test.
#!
#! Besides setting `TESTNAME`, this macro adds a new target `name` which runs the specific test alone.
macro(declaretest name)
  set(TESTNAME ${name})
  add_custom_target(${name} COMMAND ${CMAKE_CTEST_COMMAND} -V -R "^${name}$" DEPENDS ${TESTSDEPEND})
endmacro()

#! \file
#! <!--#########################################################-->
#! ### Adding all sub-directories to the test suite
#! <!--#########################################################-->
#!

add_subdirectory(run)
add_subdirectory(continue)
add_subdirectory(compile)
add_subdirectory(compilefail)
add_subdirectory(boost)
add_subdirectory(cpypyqed)
