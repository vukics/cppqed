// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "LazyDensityOperator.h"

#include <functional>



template <int RANK, typename RESULT>
using TimeDependentExpectationValue = std::function<RESULT(double t, const quantumdata::LazyDensityOperator<RANK>& psi)>;

template <int RANK, typename RESULT>
using TimeIndependentExpectationValue = std::function<RESULT(const quantumdata::LazyDensityOperator<RANK>& psi)>;



