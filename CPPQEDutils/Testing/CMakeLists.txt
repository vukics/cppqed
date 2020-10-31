# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

set(TESTING_EXECUTABLES
    HarmonicOscillator HarmonicOscillatorExact HarmonicOscillatorRandomTimes HarmonicOscillatorComplex Randomized)

foreach(S ${TESTING_EXECUTABLES})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()

enable_testing()

add_test(harmonic_oscillator_compare_ode_with_exact_solution HarmonicOscillatorComplex)
