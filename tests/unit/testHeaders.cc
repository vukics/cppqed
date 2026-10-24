// Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// Compiled twice into cppqedTests, so the framework headers that jc_qjmc uses appear in two translation units:
// any non-inline definition at namespace scope in them is a "multiple definition" at link time (#51).
#include "Composite.h"
#include "JaynesCummings.h"
#include "Mode.h"
#include "Qbit.h"
#include "Spin.h"
#include "QuantumJumpMonteCarlo.h"
#include "Pars.h"

#include <catch2/catch_test_macros.hpp>

#ifndef ODR_COPY
#define ODR_COPY "1"   // tests/CMakeLists.txt compiles this file a second time with ODR_COPY="2"
#endif

TEST_CASE("framework headers link into more than one translation unit — #51, copy " ODR_COPY, "[odr]")
{
  SUCCEED();
}
