
include_directories(${CPPQED_INCLUDE_DIRS})
message(STATUS Debug: ${CPPQED_INCLUDE_DIRS})
add_executable(ExampleModeMain ExampleModeMain.cc ExampleMode.cc ExampleModeImpl.cc ExampleInteraction.cc)
target_link_libraries(ExampleModeMain ${CPPQEDLIB})
