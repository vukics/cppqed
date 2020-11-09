# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

set(TESTING_EXECUTABLES
    EvolvedTest HarmonicOscillatorRandomTimes HarmonicOscillatorComplex MultiIndexIteratorNew Randomized)

foreach(S ${TESTING_EXECUTABLES})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()

add_test(evolved_nonlinear_ode EvolvedTest)

add_test(harmonic_oscillator_compare_ode_with_exact_solution HarmonicOscillatorComplex)

add_test(harmonic_oscillator_random_times HarmonicOscillatorRandomTimes)

add_test(multi_index_iterator MultiIndexIteratorNew)

add_test(randomized_serialization Randomized)


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

