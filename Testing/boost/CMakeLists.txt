# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

find_package(Boost ${BOOST_MINIMAL} COMPONENTS unit_test_framework)

if (NOT TARGET Boost::unit_test_framework)
  message(WARNING "No boost unit test framework found. Disabling boost tests, testing suite is not complete.")
  return()
endif()

set(BOOST_TEST_SOURCES LazyDensityOperatorFFT.cc)

add_executable(tester main.cc ${BOOST_TEST_SOURCES})
target_link_libraries(tester CPPQED::CPPQEDelements Boost::unit_test_framework)

set(BOOST_TESTS LAZY_DENSITY_OPERATOR_FFTRANSFORM)
foreach(BT ${BOOST_TESTS})
  declaretest(${BT})
  add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:tester> -t ${BT})
endforeach()
