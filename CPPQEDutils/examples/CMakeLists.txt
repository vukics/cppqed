# Copyright András Vukics 2021–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

set(EXEC_TARGETS
    HarmonicOscillator HarmonicOscillatorStdArray HarmonicOscillatorGSL HarmonicOscillatorComplexBlitz HarmonicOscillatorComplexBlitzGSL)

foreach(S ${EXEC_TARGETS})
  add_executable(${S} ${S}.cc)
  target_link_libraries(${S} CPPQED::CPPQEDutils)
endforeach()
