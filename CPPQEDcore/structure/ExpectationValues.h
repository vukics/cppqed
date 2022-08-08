// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "LazyDensityOperator.h"

#include <functional>
#include <list>
#include <valarray>
#include <variant>


namespace structure {


template <int RANK, typename RESULT>
using TimeDependentExpectationValue = std::function<RESULT(double t, const quantumdata::LazyDensityOperator<RANK>&)>;

template <int RANK, typename RESULT>
using TimeIndependentExpectationValue = std::function<RESULT(const quantumdata::LazyDensityOperator<RANK>&)>;



using EV_Array = std::valarray<double>;


/// From the size of label, it’s possible to infer the size of the array
/// pre- & postcondition: the size of EV_Array must be the same as that of label
template <int RANK>
struct ExpectationValue
{
  std::list<std::string> label;

  std::optional<std::function<void(EV_Array&)>> process{}; // by default, std::nullopt

  std::variant<TimeDependentExpectationValue<RANK,EV_Array>,TimeIndependentExpectationValue<RANK,EV_Array>> eva;

};


template <int RANK>
struct ElementExpectationValues
{
  std::string label;
  
  std::list<ExpectationValue<RANK>> expectationValues;
  
};


template <int RANK>
void streamKey(std::ostream&, const ElementExpectationValues<RANK>&, size_t);


std::ostream& stream(const EV_Array&, std::ostream& os, int precision);


/// assumption: eva[0] is expectation value, eva[1] is expectation value of the square
void calculateVariance(EV_Array& eva);


} // structure
