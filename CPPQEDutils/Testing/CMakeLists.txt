# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

set(TESTING_EXECUTABLES
    EvolvedTest HarmonicOscillatorExact HarmonicOscillatorRandomTimes HarmonicOscillatorComplex MultiIndexIteratorNew Randomized)

foreach(S ${TESTING_EXECUTABLES})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()

add_test(evolved_nonlinear_ode EvolvedTest)

add_test(harmonic_oscillator_compare_ode_with_exact_solution HarmonicOscillatorComplex)

add_test(multi_index_iterator MultiIndexIteratorNew)

add_test(randomized_serialization Randomized)


find_package(Boost ${BOOST_MINIMAL} COMPONENTS unit_test_framework)

if (TARGET Boost::unit_test_framework)
  add_executable(MultiIndexIterator MultiIndexIterator.cc)
  target_link_libraries(MultiIndexIterator CPPQED::CPPQEDutils Boost::unit_test_framework)
  add_test(multi_index_iterator_boost MultiIndexIterator)
endif()

