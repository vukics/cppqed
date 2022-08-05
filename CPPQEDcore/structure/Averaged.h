// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "LiouvillianAveragedCommon.h"

#include <list>
#include <valarray>
#include <variant>


using EV_Array = std::valarray<double>;


template <int RANK>
struct ExpectationValue
{
  std::string label; // from the size of label, it’s possible to infer the size of the array
  
  std::variant<TimeDependentExpectationValue<RANK,EV_Array>,TimeIndependentExpectationValue<RANK,EV_Array>> eva;

};


template <int RANK>
struct ElementAveraged
{
  std::string label;
  
  std::list<ExpectationValue<RANK>> expectationValues;
  
  std::optional<std::function<EV_Array&(EV_Array&)>> process;
  
};


template <int RANK>
void streamKey(std::ostream&, const ElementAveraged<RANK>&, size_t);


std::ostream& stream(const EV_Array&, std::ostream& os, int precision);

