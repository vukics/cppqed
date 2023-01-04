// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include <functional>
#include <list>
#include <valarray>


namespace structure {

using EV_Array = std::valarray<double>;

template <size_t RANK, typename RESULT=EV_Array>
using TimeDependentExpectationValue = std::function<RESULT(double t, LazyDensityOperator<RANK>)>;


template <size_t RANK, typename RESULT=EV_Array>
using TimeIndependentExpectationValue = std::function<RESULT(LazyDensityOperator<RANK>)>;


/// From the size of label, it’s possible to infer the size of the array
/// pre- & postcondition: the size of EV_Array must be the same as that of label
template <size_t RANK>
struct ExpectationValue
{
  std::list<std::string> labels; // of the ExpectationValue, like photonnumber, sigmaoperator, etc.

  std::variant<TimeDependentExpectationValue<RANK,EV_Array>,TimeIndependentExpectationValue<RANK,EV_Array>> eva;

  std::optional<std::function<void(EV_Array&)>> process{}; // by default, std::nullopt

};


template <typename T, size_t RANK>
concept expectation_values = 
  std::ranges::forward_range<T> && 
  std::is_same_v<std::ranges::range_value_t<T>,
                 ExpectationValue<RANK>>;


template <typename T>
concept key_labels = 
  std::ranges::forward_range<T> &&
  requires (std::ranges::range_value_t<T> v) {
    {v.label} -> std::same_as<std::string>;
  };
                 

template <key_labels KL>
std::ostream& streamKey(std::ostream&, const KL&);
/*{
  using namespace std;
  os<<el.label;
  for (const auto& l : el.lindblads) os<<endl<<setw(2)<<i++<<". "<<l.label;
  return os<<endl;
}*/


std::ostream& stream(const EV_Array&, std::ostream& os, int precision);


/// assumption: eva[0] is expectation value, eva[1] is expectation value of the square
void calculateVariance(EV_Array& eva);


template <size_t RANK, expectation_values<RANK> EV>
std::tuple<std::ostream&,EV_Array>
calculate_process_stream(const EV& ev, double t, const quantumdata::LazyDensityOperator<RANK>& matrix, std::ostream& os, int precision);
/*{
  Averages averages;
  if (const auto av=castAv(qs)) {
    averages.reference(av->average(t,matrix));
    av->process(averages);
    av->stream(averages,os,precision);
  }
  return {os,averages};
}*/


} // structure
