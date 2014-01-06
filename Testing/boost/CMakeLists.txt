find_package(Boost COMPONENTS unit_test_framework REQUIRED)

CPPQED_SETUP()

include(${elements_BINARY_DIR}/CPPQEDelementsConfig.cmake)
include_directories(${CPPQEDelements_INCLUDE_DIRS})

set(BOOST_TEST_SOURCES LazyDensityOperatorFFT.cc)

add_executable(tester main.cc ${BOOST_TEST_SOURCES})
target_link_libraries(tester ${CPPQED_LIBRARIES} ${CPPQEDelements_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

testdir(boost tester)

set(BOOST_TESTS LAZY_DENSITY_OPERATOR_FFTRANSFORM)
foreach(BT ${BOOST_TESTS})
  declaretest(${BT})
  add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:tester> -t ${BT})
endforeach()