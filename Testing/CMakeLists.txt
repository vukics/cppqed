file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/output")
add_test(NAME 1particle1mode
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/output"
         COMMAND 1particle1mode --dc 0 --Dt 0.05 --NDt 3 --cutoff 5 )