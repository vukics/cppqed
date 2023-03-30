// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "LazyDensityOperator.h"

#include <functional>


namespace structure {


namespace expectationvalues {

/// TODO: how to express this more conveniently (all the possibilities for a lazy_density_operator have to be listed here?)
template <typename L, size_t RANK>
concept time_dependent_functional =/*
  requires (const L& l, double t, StateVectorConstView<RANK> rho) { { l(t,rho) } -> temporal_data_point ; } &&*/
  requires (const L& l, double t, DensityOperatorConstView<RANK> rho) { { l(t,rho) } -> ::cppqedutils::temporal_data_point ; } ;

template <typename L, size_t RANK>
concept time_independent_functional =/*
  requires (const L& l, StateVectorConstView<RANK> rho) { { l(rho) } -> std::convertible_to<RESULT>; } &&*/
  requires (const L& l, DensityOperatorConstView<RANK> rho) { { l(rho) } -> ::cppqedutils::temporal_data_point ; } ;

template <typename L, size_t RANK>
concept functional = time_dependent_functional<L,RANK> || time_independent_functional<L,RANK> ;

template <size_t RANK, functional<RANK> EV>
auto calculate(const EV& ev, double t, lazy_density_operator<RANK> auto matrix)
{
  if constexpr (time_dependent_functional<EV,RANK>) return ev(t,matrix);
  else if (time_independent_functional<EV,RANK>) return ev(matrix);
}

/// pre- & postcondition: the LogTree must have the same structure as the temporal_data_point returned by expectation_value
template <typename L, size_t RANK>
concept labelled_and_without_nonlinear_postprocessing = requires (const L& l) {
  { labels(l) } -> std::convertible_to<cppqedutils::LogTree> ;
  { getFunctional(l) } -> functional<RANK> ;
};

/// Postprocessing means any operation on the expectation values that is not linear in the density operator (as the calculation of variance, for istance)
template <typename L, size_t RANK>
concept labelled_and_with_nonlinear_postprocessing = labelled_and_without_nonlinear_postprocessing<L,RANK> && (
  requires (const L& l, double t, DensityOperatorConstView<RANK> rho) { postProcess(l)( getFunctional(l)(t,rho) ); } ||
  requires (const L& l, DensityOperatorConstView<RANK> rho) { postProcess(l)( getFunctional(l)(rho) ); } ) ;

} // expectationvalues


template <typename L, size_t RANK>
concept expectation_values = expectationvalues::labelled_and_without_nonlinear_postprocessing<L,RANK> || expectationvalues::labelled_and_with_nonlinear_postprocessing<L,RANK> ;


template <size_t RANK>
void checkConsistencyWithLabels (const expectation_values<RANK> auto&);


template <size_t RANK>
auto calculateAndProcess(const expectation_values<RANK> auto& ev, double t, lazy_density_operator<RANK> auto rho) -> decltype(expectationvalues::calculate(getFunctional(ev),t,rho)) ;/*
{
  Averages averages;
  if (const auto av=castAv(qs)) {
    averages.reference(av->average(t,matrix));
    av->process(averages);
    av->stream(averages,os,precision);
  }
  return {os,averages};
}*/



/* Another definition of key_labels that mirrors that of temporal_data_points is probably an overkill
/// hana::tuple of std::string or such hana::tuples (recursive definition)
namespace expectationvalues {

template <typename T> constexpr bool kl = false;

template <> constexpr bool kl<std::string> = true;

template <hana_sequence S> constexpr bool kl<S> = !!hana::all_of(
  decltype(hana::transform(std::declval<S>(), hana::typeid_)){},
  []<class T>(T) { return kl<typename T::type>; });

} // expectationvalues

template <typename T> concept key_labels = expectationvalues::kl<T>;
*/



} // structure
