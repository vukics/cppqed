# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#! \ingroup Main
#! \file
#! \brief %CMake file for the C++QED testing suite.
#!
#! The file has the following structure:

# Please avoid any special characters in test names (+,^,.,etc.)
# The ctest command might interpret them as regular expressions


#! \file
#! <!--#########################################################-->
#! ### Configuration of the test suite
#! <!--#########################################################-->
#!
#! Find the Python interpreter and add the targets `check_fast` and `ckeck` (c.f. \ref testsuite_running).
#! Also defines helper macros.

find_package (Python 3.8 COMPONENTS Interpreter NumPy)

set ( testing_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Testing )
set ( testing_BINARY_DIR ${CMAKE_BINARY_DIR}/Testing )

# there is no official FindScipy, so we assume that NumPy implies Scipy as well
# find_package(Scipy)

if( ( TARGET Python::Interpreter ) AND ( TARGET Python::NumPy ) )
  set ( PYTHONIC_TESTDRIVER_ENABLED ON )
else()
  message(WARNING "Some tests disabled, not all Python components are available")
endif()

configure_file(testdriver.conf.in testdriver.conf)
configure_file(${CMAKE_SOURCE_DIR}/helperscripts/check_push.sh.in ${CMAKE_BINARY_DIR}/check_push.sh @ONLY)

#! \file
#! <!--#########################################################-->
#! ### Convenience custom targets
#! <!--#########################################################-->
#!

#add_custom_target(check_fast DEPENDS check_compilefail check_run check_continue check_boost testdriver.conf.in)
#add_custom_target(check DEPENDS check_fast check_compile check_cpypyqed)
#add_custom_target(check_full DEPENDS check_physics check)
#add_custom_target(check_push COMMAND ${CMAKE_BINARY_DIR}/check_push.sh)

#! \name Project options
#! @{

#! \brief Initialize a sub-directory of the test suite.
#! \ingroup Helpers
#! \param dir The sub-directory name.
#! \return `TESTSCRIPT` - set to the command by which to call the Python test driver.
#!
#! This macro adds the target `check_${dir}` and makes it depend on any target passed as argument after `dir`.
#! This target runs all tests in the current sub-directory.
#! It also sets the variable `TESTSCRIPT`, which holds the full path to the python test driver and already
#! adds some command line parameters: the configuration files and the build configuration.
#macro(testdir dir)
  #set(TESTSDEPEND ${ARGN})
  #add_custom_target(check_${dir}
                    #COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_V} ${CTEST_J}
                    #DEPENDS ${TESTSDEPEND}
  #)
  #set(TESTSCRIPT ${Python_EXECUTABLE} ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --configuration=$<CONFIGURATION>)
#endmacro()

#! \brief Declare a new test.
#! \ingroup Helpers
#! \param name Name of the test. Avoid any special characters (+,\,$ etc), which might be mis-interpreted in a regular expression.
#! \return `TESTNAME` - holds the name of the test.
#!
#! Besides setting `TESTNAME`, this macro adds a new target `name` which runs the specific test alone.
macro(declaretest name)
  set(TESTNAME ${name})
  set(TESTSCRIPT ${Python_EXECUTABLE} ${testing_SOURCE_DIR}/testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --configuration=$<CONFIGURATION>)
endmacro()

#! \file
#! <!--#########################################################-->
#! ### Adding all sub-directories to the test suite
#! <!--#########################################################-->
#!

#add_subdirectory(run)
#add_subdirectory(continue)
#add_subdirectory(compile)
add_subdirectory(compilefail)
#add_subdirectory(boost)
#add_subdirectory(cpypyqed)

if( PYTHONIC_TESTDRIVER_ENABLED )
  add_subdirectory(physics)
endif()
