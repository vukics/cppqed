// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
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


/// a label and a term
/*
 * the `label` becomes a `prefix:label` when the system becomes part of a more complex system
 */
template <size_t RANK>
using HamiltonianTerm = std::tuple<std::string,std::variant<TimeIndependentTerm<RANK>,OneTimeDependentTerm<RANK>,TwoTimeDependentTerm<RANK>>>;


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
void apply(const HA& ha, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  for (auto& term : ha)
    std::visit(cppqedutils::overload{
      [&] (const TimeIndependentTerm<RANK>& c) {c(psi,dpsidt);},
      [&] (const OneTimeDependentTerm<RANK>& c) {c(t-t0,psi,dpsidt);},
      [&] (const TwoTimeDependentTerm<RANK>& c) {c(t,psi,dpsidt,t0);},
    },std::get<1>(term));
}





/*
template <typename SYSTEM, size_t RANK>
void apply(SYSTEM&& sys, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if constexpr (std::is_convertible_v<SYSTEM,Hamiltonian<RANK>>) addContribution(std::static_cast<const >)
}
*/

} // structure

