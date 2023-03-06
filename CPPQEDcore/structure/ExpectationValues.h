// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include <functional>
#include <list>
#include <valarray>


namespace structure {


using EV_Array = std::valarray<double>;


/// TODO: how to express this more conveniently (all the possibilities for a lazy_density_operator have to be )
template <typename L, size_t RANK, typename RESULT>
concept time_dependent_expectation_value =
  requires (L&& l, double t, StateVectorConstView<RANK> rho) { { l(t,rho) } -> std::convertible_to<RESULT>; } &&
  requires (L&& l, double t, DensityOperatorConstView<RANK> rho) { { l(t,rho) } -> std::convertible_to<RESULT>; } ;

template <typename L, size_t RANK, typename RESULT>
concept time_independent_expectation_value =
  requires (L&& l, StateVectorConstView<RANK> rho) { { l(rho) } -> std::convertible_to<RESULT>; } &&
  requires (L&& l, DensityOperatorConstView<RANK> rho) { { l(rho) } -> std::convertible_to<RESULT>; } ;

template <typename L, size_t RANK, typename RESULT>
concept expectation_value = time_dependent_expectation_value<L,RANK,RESULT> || time_independent_expectation_value<L,RANK,RESULT>;


template <typename RESULT, size_t RANK, expectation_value<RANK,RESULT> EV>
RESULT calculateEV(const EV& ev, double t, lazy_density_operator<RANK> auto matrix)
{
  if constexpr (time_dependent_expectation_value<EV,RANK,RESULT>) return ev(t,matrix);
  else if (time_independent_expectation_value<EV,RANK,RESULT>) return ev(matrix);
}


template <typename T>
concept key_labels = std::ranges::forward_range<T> && std::same_as<std::ranges::range_value_t<T>,std::string>;


/// From the size of labels, it’s possible to infer the size of the array
/// pre- & postcondition: the size of EV_Array must be the same as that of label
template <typename L, size_t RANK>
concept expectation_value_element_without_process = requires (L&& l) {
 { labels(l) } -> key_labels ;
 { expectationValue(l) } -> expectation_value<RANK,EV_Array>;
};

template <typename L, size_t RANK>
concept expectation_value_element_with_process = expectation_value_element_without_process<L,RANK> && requires (L&& l) {
 { process(l) } -> std::convertible_to<std::function<void(EV_Array&)>> ;
};

template <typename L, size_t RANK>
concept expectation_value_element = expectation_value_element_without_process<L,RANK> || expectation_value_element_with_process<L,RANK> ;


template <typename L, size_t RANK>
concept expectation_values = hana_sequence<L> && !!hana::all_of(
  decltype(hana::transform(std::declval<L>(), hana::typeid_)){},
  []<class T>(T) { return expectation_value_element<typename T::type,RANK>; });


std::ostream& stream(const EV_Array&, std::ostream& os, int precision);


template <size_t RANK>
EV_Array calculateProcessStream(const expectation_values<RANK> auto& ev, double t, lazy_density_operator<RANK> auto matrix, std::ostream& os, int precision);
/*{
  Averages averages;
  if (const auto av=castAv(qs)) {
    averages.reference(av->average(t,matrix));
    av->process(averages);
    av->stream(averages,os,precision);
  }
  return {os,averages};
}*/


/*
template <size_t RANK>
struct ExpectationValue
{
  std::list<std::string> labels; // of the ExpectationValue, like photonnumber, sigmaoperator, etc.

  std::variant<TimeDependentExpectationValue<RANK,EV_Array>,TimeIndependentExpectationValue<RANK,EV_Array>> eva;

  std::optional<> process{}; // by default, std::nullopt

};*/

                 

template <key_labels KL>
std::ostream& streamKey(std::ostream&, const KL&);
/*{
  using namespace std;
  os<<el.label;
  for (const auto& l : el.lindblads) os<<endl<<setw(2)<<i++<<". "<<l.label;
  return os<<endl;
}*/




/// assumption: eva[0] is expectation value, eva[1] is expectation value of the square
void calculateVariance(EV_Array& eva);


} // structure
