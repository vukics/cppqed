// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExpectationValues.h"


namespace structure {

  
using Rates=EV_Array;


template <size_t RANK>
using TimeDependentJump = std::function<void(double t, StateVectorView<RANK> psi)>;

template <size_t RANK>
using TimeIndependentJump = std::function<void(StateVectorView<RANK> psi)>;


template <size_t RANK>
using TimeDependentRate = TimeDependentExpectationValue<RANK,double>;

template <size_t RANK>
using TimeIndependentRate = TimeIndependentExpectationValue<RANK,double>;


template <size_t RANK>
using TimeDependentSuperoperator = std::function<void(double t, DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK>& drhodt)>;

template <size_t RANK>
using TimeIndependentSuperoperator = std::function<void(DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK> drhodt)>;


template <size_t RANK>
using Rate = std::variant<TimeDependentRate<RANK>,TimeIndependentRate<RANK>>;

template <size_t RANK>
using Jump = std::variant<TimeDependentJump<RANK>,TimeIndependentJump<RANK>>;

template <size_t RANK>
using Superoperator = std::variant<TimeDependentSuperoperator<RANK>,TimeIndependentSuperoperator<RANK>>;


template <size_t RANK>
struct Lindblad
{
  std::string label;
  Jump<RANK> jump;
  Rate<RANK> rate; // could be made std::optional<Rate<RANK>> rate; if it’s null, then the rate is calculated from the jump functional
  Superoperator<RANK> superoperator; // std::optional<Superoperator<RANK>> superoperator; if it’s null, then the superoperator is calculated from the jump functional
};


template <typename T, size_t RANK>
concept liouvillian = 
  std::ranges::forward_range<T> && 
  std::is_same_v<std::ranges::range_value_t<T>,
                 Lindblad<RANK>>;


template <size_t RANK, liouvillian<RANK> LI>
EV_Array calculateRates(const LI& li, double t, LazyDensityOperator<RANK> matrix);
// {
//   for (auto& term : ha)
//     std::visit(cppqedutils::overload{
//       [&] (const TimeIndependentTerm<RANK>& c) {c(psi,dpsidt);},
//       [&] (const OneTimeDependentTerm<RANK>& c) {c(t-t0,psi,dpsidt);},
//       [&] (const TwoTimeDependentTerm<RANK>& c) {c(t,psi,dpsidt,t0);},
//     },std::get<1>(term));
// }

                 

/// alternatively, we could use this and the next function as default values for rate and superoperator in Lindblads, but this seems a bit over-the-top
template <size_t RANK>
auto rateFromJump(double t, const quantumdata::StateVector<RANK>& psi, Jump<RANK> jump)/*
{
  quantumdata::StateVector<RANK> psiTemp(psi);
  if (!jump.index()) std::get<0>(jump)(t,psiTemp.getArray());
  else std::get<1>(jump)(psiTemp.getArray());
  return std::tuple{cppqedutils::sqr(psiTemp.norm()),
                    psiTemp // careful! psiTemp is not normalized
                   };
}*/;


template <size_t RANK>
void superoperatorFromJump(double t, DensityOperatorConstView<RANK> rho, DensityOperatorView<RANK> drhodt, Jump<RANK> jump)/*
{
  DensityOperatorConstView<RANK> rhotemp(rho.copy());

  auto unaryIteration=[&] () {
    for (auto& psi : fullRange<blitzplusplus::vfmsi::Left>(rhotemp)) {
      if (!jump.index()) std::get<0>(jump)(t,psi);
      else std::get<1>(jump)(psi);
    }
  };
  
  unaryIteration();
  blitzplusplus::hermitianConjugateSelf(rhotemp);
  unaryIteration();

  drhodt+=rhotemp;
}*/;



} // structure
