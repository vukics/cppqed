// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include <functional>


namespace structure {

using namespace quantumdata;

namespace expectation_values_ns {

static const struct NoOp
{
  std::string label{"noOp"};
  
  auto operator() (auto) const {return hana::make_tuple();}

} noOp;


template <typename L, size_t RANK>
concept time_dependent_functional = requires (const L& l, LDO<StateVector,RANK> rho, double t)
{
  { l(t,rho) } -> temporal_data_point ;
};


template <typename L, size_t RANK>
concept time_independent_functional = requires (const L& l, LDO<StateVector,RANK> rho)
{ 
  { l(rho) } -> temporal_data_point ;
};


// Alternatively:
/* requires (const L& l) { requires (
requires (StateVectorConstView<RANK> rho) { { l(rho) } -> temporal_data_point ; } &&
requires (DensityOperatorConstView<RANK> rho) { { l(rho) } -> temporal_data_point ; } ) ; } ;*/

template <typename L, size_t RANK>
concept functional = time_dependent_functional<L,RANK> || time_independent_functional<L,RANK> ;


} // expectation_values_ns


/// pre- & postcondition: the LogTree must have the same structure as the temporal_data_point returned by expectation_value
template <typename L, size_t RANK>
concept expectation_values = /* labelled<L> && */ expectation_values_ns::functional<L,RANK> ;


template <size_t RANK>
void checkConsistencyWithLabels (const expectation_values<RANK> auto&);


template <size_t RANK, expectation_values<RANK> EV>
auto calculateExpectationValues(const EV& ev, double t, lazy_density_operator<RANK> auto matrix)
{
  if constexpr (expectation_values_ns::time_dependent_functional<EV,RANK>) return ev(t,matrix);
  else return ev(matrix);
}


namespace expectation_values_ns {

template <
  auto retainedAxes,
  size_t RANK,
  functional<std::size(retainedAxes)> EV>
auto broadcast(const EV& ev, double t, lazy_density_operator<RANK> auto matrix, const std::vector<size_t>& offsets)
{
  return partialTrace( matrix, offsets, [&] (auto psiElem) {return calculate(ev,t,psiElem); } );
}

} // expectation_values_ns


/* Another definition of key_labels that mirrors that of temporal_data_points is probably an overkill
/// hana::tuple of std::string or such hana::tuples (recursive definition)
namespace expectation_values_ns {

template <typename T> constexpr bool kl = false;

template <> constexpr bool kl<std::string> = true;

template <hana_sequence S> constexpr bool kl<S> = !!hana::all_of(
  decltype(hana::transform(std::declval<S>(), hana::typeid_)){},
  []<class T>(T) { return kl<typename T::type>; });

} // expectation_values_ns

template <typename T> concept key_labels = expectation_values_ns::kl<T>;
*/



} // structure
