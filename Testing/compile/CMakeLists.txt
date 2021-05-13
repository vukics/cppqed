add_library(SliceIteratorQuantumData OBJECT SliceIteratorQuantumData.cc)
target_link_libraries(SliceIteratorQuantumData CPPQED::CPPQEDcore_quantumdata)
add_test( NAME test_SliceIteratorQuantumData COMMAND cmake --build ${CMAKE_CURRENT_BINARY_DIR} --target SliceIteratorQuantumData )
set_tests_properties(test_SliceIteratorQuantumData PROPERTIES LABELS "Fast")

