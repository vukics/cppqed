# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_executable(ExampleModeMain ExampleModeMain.cc ExampleMode.cc ExampleModeImpl.cc ExampleInteraction.cc)
target_link_libraries(ExampleModeMain CPPQED::CPPQEDcore)

#add_executable(HarmonicOscillatorComplex HarmonicOscillatorComplex.cc)
#target_link_libraries(HarmonicOscillatorComplex CPPQED::CPPQEDcore)

add_executable(CoherentElementTest CoherentElementTest.cc)
target_link_libraries(CoherentElementTest CPPQED::CPPQEDcore)

add_library(ExampleObjects OBJECT cpcPaperExamples.cc)
target_link_libraries(ExampleObjects CPPQED::CPPQEDcore)
