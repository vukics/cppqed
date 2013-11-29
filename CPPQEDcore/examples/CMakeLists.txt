include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
CPPQED_SETUP()
message(STATUS "Debug: ${CPPQED_INCLUDE_DIRS}")
add_executable(ExampleModeMain ExampleModeMain.cc ExampleMode.cc ExampleModeImpl.cc ExampleInteraction.cc)
target_link_libraries(ExampleModeMain ${CPPQED_LIBRARIES})

add_executable(CoherentElementTest CoherentElementTest.cc)
target_link_libraries(CoherentElementTest ${CPPQED_LIBRARIES})
