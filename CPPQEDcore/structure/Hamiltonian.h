// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"

#include "Overloaded.h"

#include <variant>


namespace structure {


template <size_t RANK>
using TimeIndependentTerm = std::function<void(StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt)>;


template <size_t RANK>
using OneTimeDependentTerm = std::function<void(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt)>;


template <size_t RANK>
using TwoTimeDependentTerm = std::function<void(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)>;


template <size_t RANK>
using OneTimeDependentPropagator = std::function<void(double, StateVectorView<RANK>)>;


template <size_t RANK>
using TwoTimeDependentPropagator = std::function<void(double, StateVectorView<RANK>, double)>;


/// a label and a term
/*
 * the `label` becomes a `prefix:label` when the system becomes part of a more complex system
 */
template <size_t RANK>
using HamiltonianTerm = std::tuple<std::string,std::variant<TimeIndependentTerm<RANK>,OneTimeDependentTerm<RANK>,TwoTimeDependentTerm<RANK>,OneTimeDependentPropagator<RANK>,TwoTimeDependentPropagator<RANK>>>;


// e.g. `std::list<HamiltonianTerm<RANK>>`, but transformed views might be necessary
template <typename T, size_t RANK>
concept hamiltonian = 
  std::ranges::forward_range<T> && 
  std::is_same_v<std::ranges::range_value_t<T>,
                 HamiltonianTerm<RANK>>;


/// applying a Hamiltonian is by default interpreted as |dpsidt>+=H|psi>
/**
 * However, if indexing is done carefully, `psi` and `dpsidt` can refer to the same underlying data
 */
template <size_t RANK, hamiltonian<RANK> HA>
void applyHamiltonian(const HA& ha, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  for (auto& term : ha)
    std::visit(cppqedutils::overload{
      [&] (const TimeIndependentTerm<RANK>& c) {c(psi,dpsidt);},
      [&] (const OneTimeDependentTerm<RANK>& c) {c(t-t0,psi,dpsidt);},
      [&] (const TwoTimeDependentTerm<RANK>& c) {c(t,psi,dpsidt,t0);},
    },std::get<1>(term));
}


/// applying a Propagator is interpreted as replacing |psi> with U|psi>
template <size_t RANK, hamiltonian<RANK> HA>
void applyPropagator(const HA& ha, double t, StateVectorView<RANK> psi, double t0)
{
  for (auto& term : ha)
    std::visit(cppqedutils::overload{
      [&] (const OneTimeDependentPropagator<RANK>& c) {c(t-t0,psi);},
      [&] (const TwoTimeDependentPropagator<RANK>& c) {c(t,psi,t0);},
    },std::get<1>(term));
}



/// A unary propagator that assumes that the operator that transforms between the pictures is diagonal
template<bool IS_TWO_TIME>
struct UnaryDiagonalPropagator
{
  using Diagonal = std::valarray<dcomp>;

  using UpdateFunctional = std::conditional_t<IS_TWO_TIME,std::function<void(double,double,Diagonal&)>,std::function<void(double,Diagonal&)>>;
  
  UnaryDiagonalPropagator(size_t dim, UpdateFunctional updateDiagonal) : diagonal{dim}, updateDiagonal_{updateDiagonal} {}
  
  void operator()(double t, StateVectorView<1> psi, double t0) const requires (IS_TWO_TIME) {if (t!=t_ || t0!=t0_) {updateDiagonal_(t_=t,t0_=t0,diagonal);} psi*=diagonal;}

  void operator()(double t, StateVectorView<1> psi) const requires (!IS_TWO_TIME) {if (t!=t_) {updateDiagonal_(t_=t,diagonal);} psi*=diagonal;}

  mutable Diagonal diagonal;

private:
  mutable double t_=0, t0_=0;
  
  const UpdateFunctional updateDiagonal_;

};



} // structure

