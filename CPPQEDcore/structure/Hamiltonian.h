// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"


namespace structure {

using namespace quantumdata;

/// TODO: this could be fused with ode::system, but let’s not go into that for the moment :)
template <typename H, typename StateIn, typename StateOut>
concept ode_time_independent_derivative = requires(const H& h, StateIn psi, StateOut dpsidt) { h(psi,dpsidt); };

template <typename H, typename StateIn, typename StateOut>
concept ode_time_dependent_derivative = requires(const H& h, double t, StateIn psi, StateOut dpsidt) { h(t,psi,dpsidt); };


namespace hamiltonian_ns {

static const struct NoOp {LogTree label{"noOp"};} noOp;

template <typename H, size_t RANK>
concept time_independent_functional = ode_time_independent_derivative<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept one_time_dependent_functional = ode_time_dependent_derivative<H,StateVectorConstView<RANK>,StateVectorView<RANK>>;

template <typename H, size_t RANK>
concept two_time_dependent_functional = requires (const H& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) { h(t,psi,dpsidt,t0); };

template <typename H, size_t RANK>
concept functional = time_independent_functional<H,RANK> || one_time_dependent_functional<H,RANK> || two_time_dependent_functional<H,RANK> || std::same_as<std::decay_t<H>,NoOp>;


} // hamiltonian_ns


/// applying a Hamiltonian term is by default interpreted as |dpsidt>+=H|psi>/(i*hbar)
/** However, if indexing is done carefully, `psi` and `dpsidt` can refer to the same underlying data */
template <size_t RANK, hamiltonian_ns::functional<RANK> T>
void applyHamiltonian(const T& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if      constexpr (hamiltonian_ns::  time_independent_functional<T,RANK>) h(psi,dpsidt);
  else if constexpr (hamiltonian_ns::one_time_dependent_functional<T,RANK>) h(t-t0,psi,dpsidt);
  else if constexpr (hamiltonian_ns::two_time_dependent_functional<T,RANK>) h(t,psi,dpsidt,t0);
}


/// composition of two hamiltonian functionals, meant for interaction-type elements
/**
 * For simplicity, it always returns a hamiltonian_ns::two_time_dependent_functional
 * \note Here, we cannot really use the operator syntax as in MultiDiagonal as this function will not be found by ADL
 */
template <size_t RANK1, size_t RANK2>
auto compose(const hamiltonian_ns::functional<RANK1> auto& h1, const hamiltonian_ns::functional<RANK2> auto& h2);




template <typename H, size_t RANK>
concept hamiltonian = labelled<H> && hamiltonian_ns::functional<H,RANK>;



template <size_t RANK, hamiltonian_ns::functional<RANK> F>
struct HamiltonianElement
{
  std::string label;

  F functional;
  
  void operator()(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) const
  {
    applyHamiltonian(functional,t,psi,dpsidt,t0);
  };

};


/// TODO: how to make a deduction guide that would deduce RANK from the StateVector argument of the functional?
template <size_t RANK, hamiltonian_ns::functional<RANK> F>
auto makeHamiltonianElement(std::string label, F&& functional)
{
  return HamiltonianElement<RANK,F>{.label{label},.functional{std::forward<F>(functional)}};
}


// itself a hamiltonian
/**
 * \note It would be easy to do a runtime collection as well. For that, all the functionals should be wrapped in a `std::function`, and stored in a runtime container
 */
template <size_t RANK, hamiltonian<RANK>... H>
struct HamiltonianCollection
{
  hana::tuple<H...> collection;

  void operator()(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0) const
  {
    hana::for_each(collection, [&] (const auto& h) {applyHamiltonian(h,t,psi,dpsidt,t0);});
  }

  friend LogTree label(const HamiltonianCollection& hc) {
    LogTree res;
    hana::for_each(hc.collection,[&] (const auto& h) {res.push_back(getLabel(h));});
    return res;
  }

};


template <size_t RANK, hamiltonian<RANK>... H>
auto makeHamiltonianCollection(H&&... h) { return HamiltonianCollection<RANK,H...>{.collection{std::forward<H>(h)...}}; }



/*
template <typename H, size_t RANK>
concept hamiltonian = hana_sequence<H> && !!hana::all_of(
  decltype(hana::transform(std::declval<H>(), hana::typeid_)){},
  []<class T>(T) { return hamiltonian_element<typename T::type,RANK>; });

*/

namespace hamiltonian_ns {

/// Maybe this can be just an overload of applyHamiltonian ?
// TODO: std::ranges::views::zip to be applied here, but for some reason, sliceRange is not compatible with zipping – it does work with sliceRangeSimple, though.
//for (auto [psi,dpsidt] : std::views::zip(sliceRange<retainedAxes>(psi,offsets),
//                                         sliceRange<retainedAxes>(drhodt,offsets) ) )
//  ::structure::applyHamiltonian(h,t,psi,dpsidt,t0) ;
template <
  auto retainedAxes,
  size_t RANK,
  functional<std::size(retainedAxes)> T >
void broadcast(const T& h, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0, const std::vector<size_t>& offsets)
{
  for ( auto&& [psi,dpsidt] : std::views::zip( sliceRange<retainedAxes>(psi,offsets), sliceRange<retainedAxes>(dpsidt,offsets) ) )
    applyHamiltonian(h,t,psi,dpsidt,t0);
}


} // hamiltonian_ns


template <typename StateIn, typename StateOut>
using ODE_derivativeTimeIndependentFunctional = std::function<void(StateIn psi, StateOut dpsidt)> ;

template <typename StateIn, typename StateOut>
using ODE_derivativeTimeDependentFunctional = std::function<void(double t, StateIn psi, StateOut dpsidt)> ;


template <size_t RANK> using TimeIndependentTerm = ODE_derivativeTimeIndependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;
template <size_t RANK> using TimeDependentTerm = ODE_derivativeTimeDependentFunctional<StateVectorConstView<RANK>,StateVectorView<RANK>>;


} // structure

