# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

find_package(Boost ${BOOST_MINIMAL} COMPONENTS unit_test_framework)

if (NOT TARGET Boost::unit_test_framework)
  message(WARNING "No boost unit test framework found. Disabling boost tests, testsuite is not complete.")
  return()
endif()

set(BOOST_TEST_SOURCES LazyDensityOperatorFFT.cc EntanglementMeasures.cc)

add_executable(boostTester main.cc ${BOOST_TEST_SOURCES})
target_link_libraries(boostTester CPPQED::CPPQEDelements Boost::unit_test_framework)

set(BOOST_TESTS LAZY_DENSITY_OPERATOR_FFTRANSFORM NEGATIVITY_OF_PARTIAL_TRANSPOSE)
foreach(BT ${BOOST_TESTS})
  declaretest(${BT})
  add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:boostTester> -t ${BT})
endforeach()
