# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
# RATIONALE:
# * any test that runs and checks for a single case, just uses plain ctest here
# * those that use more cases rely on boost::unit_test_framework

macro(add_test_utils testName testTarget)
  add_test(${testName} ${testTarget} ${ARGN})
  set_tests_properties(${testName} PROPERTIES LABELS "Fast;Utils")
endmacro()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

set(TESTING_EXECUTABLES
    EvolvedTest EvolvedTestGSL MultiIndexIteratorNew)

if (TARGET Eigen3::Eigen)
  set(TESTING_EXECUTABLES ${TESTING_EXECUTABLES} HarmonicOscillatorRandomTimes HarmonicOscillatorComplex HarmonicOscillatorComplexGSL)
endif()
    
foreach(S ${TESTING_EXECUTABLES})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()

add_test_utils(evolved_nonlinear_ode_boost EvolvedTest)

add_test_utils(evolved_nonlinear_ode_gsl EvolvedTestGSL)


if (TARGET Eigen3::Eigen)
  add_test_utils(harmonic_oscillator_compare_ode_with_exact_solution HarmonicOscillatorComplex)
  add_test_utils(harmonic_oscillator_compare_ode_with_exact_solution_gsl HarmonicOscillatorComplexGSL)
  add_test_utils(harmonic_oscillator_random_times HarmonicOscillatorRandomTimes)
endif()

add_test_utils(multi_index_iterator MultiIndexIteratorNew)

find_package(Boost ${BOOST_MINIMAL} COMPONENTS unit_test_framework)

if (TARGET Boost::unit_test_framework)
  set(BOOST_TESTING_EXECUTABLES
      BlitzArray BlitzExtensions BlitzTinyExtensions CWR MultiIndexIterator Random SliceIterator)
  foreach(S ${BOOST_TESTING_EXECUTABLES})
    add_executable(${S} ${S}.cc)
    target_link_libraries(${S} CPPQED::CPPQEDutils Boost::unit_test_framework)
  endforeach()
  
  add_test_utils(blitz_array BlitzArray)
  
  add_test_utils(blitz_extensions BlitzExtensions)
  
  add_test_utils(blitz_tiny_extensions BlitzTinyExtensions)
  
  add_test_utils(combination_with_repetition CWR)
  
  add_test_utils(multi_index_iterator_boost MultiIndexIterator)

  add_test_utils(random_engine_serialization ${CMAKE_CURRENT_BINARY_DIR}/Random)
  set_tests_properties(random_engine_serialization PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

  add_test_utils(slice_iterator SliceIterator)
endif()

set(COMPILE_TESTING_NAMES compile_BlitzArraySliceIteratorTMP compile_TMP_Tools)

foreach(S ${COMPILE_TESTING_NAMES})
  add_library(${S} OBJECT ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
  add_test( NAME test_${S} COMMAND cmake --build ${CMAKE_CURRENT_BINARY_DIR} --target ${S} )
  set_tests_properties(test_${S} PROPERTIES LABELS "Fast;Utils")
  # add_test( NAME ${S} COMMAND ${PYTHON_EXECUTABLE} testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --test=${TESTNAME} --testclass=CompileTarget --script=TMP_ToolsCompile)
endforeach()

macro(add_compilefail_test S errorText)
  add_library(${S} OBJECT ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils ${ARGN})
  set_target_properties(${S} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  add_test( NAME test_${S} COMMAND cmake --build ${CMAKE_CURRENT_BINARY_DIR} --target ${S} )
  set_tests_properties(test_${S} PROPERTIES PASS_REGULAR_EXPRESSION ${errorText} LABELS "Fast")
endmacro()


macro(add_compilefail_test_utils S errorText)
 add_compilefail_test(${S} ${errorText})
 set_tests_properties(test_${S} PROPERTIES LABELS "Utils")
endmacro()


add_compilefail_test_utils(fail_BlitzArraySliceIterator1 "cppqedutils::SliceIterator inconsistent vector")

add_compilefail_test_utils(fail_BlitzArraySliceIterator2 "cppqedutils::SliceIterator vector out of range")

add_compilefail_test_utils(fail_BlitzTinyExtensions "candidate template ignored: requirement '7 % 2 == 0' was not satisfied")

add_compilefail_test_utils(fail_TMP_Tools1 "Negative element in nonnegative vector")

# add_compilefail_test_utils(fail_TMP_Tools2 "Argument not even") # obsolate

add_compilefail_test_utils(fail_TMP_Tools3 "pair_c with equal elements")


add_custom_target(compileTests DEPENDS ${TESTING_EXECUTABLES} ${BOOST_TESTING_EXECUTABLES})


# find_package (Python 3 COMPONENTS Interpreter NumPy)
# find_package(Scipy)
