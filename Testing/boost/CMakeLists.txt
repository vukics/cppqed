# Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

find_package(Boost QUIET COMPONENTS unit_test_framework)
add_feature_info("Boost unit test framework"
                 Boost_UNIT_TEST_FRAMEWORK_FOUND
                 "Part of the test suite." )

if(NOT Boost_UNIT_TEST_FRAMEWORK_FOUND)
  message(WARNING "No boost unit test framework found. Disabling boost tests, testing suite is not complete.")
  return()
endif()

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