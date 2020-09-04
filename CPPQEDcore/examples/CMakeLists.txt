include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
CPPQED_SETUP()
add_executable(ExampleModeMain ExampleModeMain.cc ExampleMode.cc ExampleModeImpl.cc ExampleInteraction.cc)
target_link_libraries(ExampleModeMain ${CPPQED_LIBRARIES})

add_executable(HarmonicOscillatorComplex HarmonicOscillatorComplex.cc)
target_link_libraries(HarmonicOscillatorComplex ${CPPQED_LIBRARIES})

add_executable(CoherentElementTest CoherentElementTest.cc)
target_link_libraries(CoherentElementTest ${CPPQED_LIBRARIES})

add_library(ExampleObjects OBJECT cpcPaperExamples.cc)

foreach(target ExampleModeMain HarmonicOscillatorComplex CoherentElementTest ExampleObjects)
  if(BUNDLED_BLITZ)
    add_dependencies(${target} blitz)
  endif()
endforeach()

add_custom_target(core_examples_main_targets)
add_dependencies(core_examples_main_targets ExampleModeMain HarmonicOscillatorComplex CoherentElementTest)
