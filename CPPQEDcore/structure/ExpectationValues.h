// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include <functional>
#include <list>
#include <valarray>


namespace structure {

using EV_Array = std::valarray<double>;

template <typename L, size_t RANK, typename LDO, typename RESULT>
concept time_dependent_expectation_value = lazy_density_operator<LDO,RANK> && requires (L&& l, double t, LDO rho) { { l(t,rho) } -> std::convertible_to<RESULT>; } ;

template <typename L, size_t RANK, typename LDO, typename RESULT>
concept time_independent_expectation_value = lazy_density_operator<LDO,RANK> && requires (L&& l, LDO rho) { { l(rho) } -> std::convertible_to<RESULT>; } ;

template <typename L, size_t RANK, typename LDO, typename RESULT>
concept expectation_value = time_dependent_expectation_value<L,RANK,LDO,RESULT> || time_independent_expectation_value<L,RANK,LDO,RESULT>;


template <size_t RANK, lazy_density_operator<RANK> LDO, typename RESULT, expectation_value<RANK,LDO,RESULT> EV>
EV_Array calculateEV(const EV& ev, double t, LDO matrix)
{
  if constexpr (time_dependent_expectation_value<EV,RANK,LDO,RESULT>) return ev(t,matrix);
  else if (time_independent_expectation_value<EV,RANK,LDO,RESULT>) return ev(matrix);
}


template <typename T>
concept key_labels = std::ranges::forward_range<T> && std::same_as<std::ranges::range_value_t<T>,std::string>;


/// From the size of labels, it’s possible to infer the size of the array
/// pre- & postcondition: the size of EV_Array must be the same as that of label
template <typename L, size_t RANK, typename LDO>
concept expectation_value_element_without_process = requires (L&& l) {
 { labels(l) } -> key_labels ;
 { expectationValue(l) } -> expectation_value<RANK,LDO,EV_Array>;
};

template <typename L, size_t RANK, typename LDO>
concept expectation_value_element_with_process = expectation_value_element_without_process<L,RANK,LDO> && requires (L&& l) { 
 { process(l) } -> std::convertible_to<std::function<void(EV_Array&)>> ;
};

template <typename L, size_t RANK, typename LDO>
concept expectation_value_element = expectation_value_element_without_process<L,RANK,LDO> || expectation_value_element_with_process<L,RANK,LDO> ;


template <typename L, size_t RANK, typename LDO>
concept expectation_values = hana_sequence<L> && !!hana::all_of(
  decltype(hana::transform(std::declval<L>(), hana::typeid_)){},
  []<class T>(T) { return expectation_value_element<typename T::type,RANK,LDO>; });


std::ostream& stream(const EV_Array&, std::ostream& os, int precision);



template <size_t RANK, lazy_density_operator<RANK> LDO>
EV_Array calculateProcessStream(const expectation_values<RANK,LDO> auto& ev, double t, LDO matrix, std::ostream& os, int precision);
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
