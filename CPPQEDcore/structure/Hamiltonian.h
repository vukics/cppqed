// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "StateVector.h"

#include "Overoaded.h"

#include <variant>


namespace structure {


template <size_t RANK>
using TimeIndependentContribution = std::function<void(StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt)>;


template <size_t RANK>
using OneTimeDependentContribution = std::function<void(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt)>;


template <size_t RANK>
using TwoTimeDependentContribution = std::function<void(double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)>;


template <size_t RANK>
using HamiltonianContribution = std::variant<TimeIndependentContribution<RANK>,OneTimeDependentContribution<RANK>,TwoTimeDependentContribution<RANK>>;


template <size_t RANK>
using Hamiltonian = std::list<HamiltonianContribution<RANK>>;


template <typename H, size_t RANK>
concept HamiltonianSystem = std::convertible_to<H,HamiltonianContribution<RANK>>;


template <size_t RANK>// requires Hamiltonian<H,RANK>
void addContribution(const Hamiltonian<RANK>& hamiltonian, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  for (auto& contribution : hamiltonian)
    std::visit(cppqedutils::overload{
      [&] (const TimeIndependentContribution<RANK>& c) {c(psi,dpsidt);},
      [&] (const OneTimeDependentContribution<RANK>& c) {c(t-t0,psi,dpsidt);},
      [&] (const TwoTimeDependentContribution<RANK>& c) {c(t,psi,dpsidt,t0);},
    },contribution);
}





/*
template <typename SYSTEM, size_t RANK>
void addContribution(SYSTEM&& sys, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if constexpr (std::is_convertible_v<SYSTEM,Hamiltonian<RANK>>) addContribution(std::static_cast<const >)
}
*/

} // structure

