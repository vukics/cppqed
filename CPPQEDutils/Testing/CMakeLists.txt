# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
# RATIONALE:
# * any test that runs and checks for a single case, just uses plain ctest here
# * those that use more cases rely on boost::unit_test_framework

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

configure_file(testdriver.conf.in testdriver.conf)

set(TESTING_EXECUTABLES
    EvolvedTest MultiIndexIteratorNew Randomized)

if (TARGET Eigen3::Eigen)
  set(TESTING_EXECUTABLES ${TESTING_EXECUTABLES} HarmonicOscillatorRandomTimes HarmonicOscillatorComplex)
endif()
    
foreach(S ${TESTING_EXECUTABLES})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()

add_test(evolved_nonlinear_ode EvolvedTest)

if (TARGET Eigen3::Eigen)
  add_test(harmonic_oscillator_compare_ode_with_exact_solution HarmonicOscillatorComplex)
  add_test(harmonic_oscillator_random_times HarmonicOscillatorRandomTimes)
endif()

add_test(multi_index_iterator MultiIndexIteratorNew)

add_test(randomized_serialization ${CMAKE_CURRENT_BINARY_DIR}/Randomized)
set_tests_properties(randomized_serialization PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

find_package(Boost ${BOOST_MINIMAL} COMPONENTS unit_test_framework)

if (TARGET Boost::unit_test_framework)
  set(BOOST_TESTING_EXECUTABLES
      BlitzArray BlitzExtensions BlitzTinyExtensions CWR MultiIndexIterator SliceIterator)
  foreach(S ${BOOST_TESTING_EXECUTABLES})
    add_executable(${S} ${S}.cc)
    target_link_libraries(${S} CPPQED::CPPQEDutils Boost::unit_test_framework)
  endforeach()
  
  add_test(blitz_array BlitzArray)
  
  add_test(blitz_extensions BlitzExtensions)
  
  add_test(blitz_tiny_extensions BlitzTinyExtensions)
  
  add_test(combination_with_repetition CWR)
  
  add_test(multi_index_iterator_boost MultiIndexIterator)
  
  add_test(slice_iterator SliceIterator)
endif()

set(COMPILE_TESTING_NAMES compile_BlitzArraySliceIteratorTMP compile_TMP_Tools)

foreach(S ${COMPILE_TESTING_NAMES})
  add_library(${S} OBJECT ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
  add_test( NAME test_${S} COMMAND cmake --build ${CMAKE_CURRENT_BINARY_DIR} --target ${S} )
  # add_test( NAME ${S} COMMAND ${PYTHON_EXECUTABLE} testdriver.py ${testing_BINARY_DIR}/testdriver.conf ${CMAKE_CURRENT_LIST_DIR}/testdriver.conf --test=${TESTNAME} --testclass=CompileTarget --script=TMP_ToolsCompile)
endforeach()

macro(add_compilefail_test S errorText)
  add_library(${S} OBJECT ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
  set_target_properties(${S} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  add_test( NAME test_${S} COMMAND cmake --build ${CMAKE_CURRENT_BINARY_DIR} --target ${S} )
  set_tests_properties(test_${S} PROPERTIES PASS_REGULAR_EXPRESSION ${errorText})
endmacro()

add_compilefail_test(fail_BlitzArraySliceIterator1 "cpputils::SliceIterator inconsistent vector")

add_compilefail_test(fail_BlitzArraySliceIterator2 "cpputils::SliceIterator vector out of range")

# ??? fail_BlitzTinyExtensions compiles without error ???
# add_compilefail_test(fail_BlitzTinyExtensions "Argument not even")

add_compilefail_test(fail_TMP_Tools1 "Negative element in nonnegative vector")

add_compilefail_test(fail_TMP_Tools2 "Argument not even")

add_compilefail_test(fail_TMP_Tools3 "pair_c with equal elements")


add_custom_target(compileTests DEPENDS ${TESTING_EXECUTABLES} ${BOOST_TESTING_EXECUTABLES})

# find_package (Python 3 COMPONENTS Interpreter NumPy)

# find_package(Scipy)
