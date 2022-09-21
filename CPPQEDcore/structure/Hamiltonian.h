// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "StateVector.h"

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


template <typename H, size_t RANK>
concept Hamiltonian = std::ranges::input_range<H> && std::convertible_to<std::ranges::range_value_t<H>,HamiltonianContribution<RANK>>;


template <typename H, size_t RANK> requires Hamiltonian<H,RANK>
void addContribution(const H& hamiltonian, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  for (auto&& convertible : hamiltonian) {
    HamiltonianContribution<RANK> contribution{convertible};
    switch (contribution.index()) {
      case 0 : std::get<0>(contribution)(psi,dpsidt); break;
      case 1 : std::get<1>(contribution)(t-t0,psi,dpsidt); break;
      case 2 : std::get<2>(contribution)(t,psi,dpsidt,t0); break;
    }
  }
}


/*
template <typename SYSTEM, size_t RANK>
void addContribution(SYSTEM&& sys, double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt, double t0)
{
  if constexpr (std::is_convertible_v<SYSTEM,Hamiltonian<RANK>>) addContribution(std::static_cast<const >)
}
*/

} // structure

