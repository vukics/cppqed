// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExpectationValues.h"
#include "Hamiltonian.h"


namespace structure {

  
using Rates=EV_Array;


template <typename L, size_t RANK>
concept time_dependent_jump = one_time_dependent_propagator<L,RANK>;

template <typename L, size_t RANK>
concept time_independent_jump = requires(L&& l, StateVectorView<RANK> psi) { l(psi); };

template <typename L, size_t RANK>
concept jump = time_dependent_jump<L,RANK> || time_independent_jump<L,RANK>;


template <typename L, size_t RANK>
concept time_dependent_rate = time_dependent_expectation_value<L,RANK,double>;

template <typename L, size_t RANK>
concept time_independent_rate = time_independent_expectation_value<L,RANK,double>;

template <typename L, size_t RANK>
concept rate = time_dependent_rate<L,RANK> || time_independent_rate<L,RANK>;


template <typename L, size_t RANK>
concept time_dependent_superoperator = ode_derivative_time_dependent_contribution<L,DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;

template <typename L, size_t RANK>
concept time_independent_superoperator = ode_derivative_time_independent_contribution<L,DensityOperatorConstView<RANK>,DensityOperatorView<RANK>>;

template <typename L, size_t RANK>
concept superoperator = time_dependent_superoperator<L,RANK> || time_independent_superoperator<L,RANK>;


template <typename L, size_t RANK>
concept lindblad_with_jump = labelled<L> && requires (L&& l) {
  { jumpFunctional(l) } -> jump<RANK> ;
};

template <typename L, size_t RANK>
concept lindblad_with_rate = lindblad<L> && requires (L&& l) {
  { rateFunctional(l) } -> rate<RANK> ;
};

template <typename L, size_t RANK>
concept lindblad_with_superoperator = lindblad<L> && requires (L&& l) {
  { superoperatorFunctional(l) } -> superoperator<RANK> ;
};

template <typename L, size_t RANK>
concept lindblad = lindblad_with_superoperator<L,RANK> || lindblad_with_rate<L,RANK> || lindblad_with_jump<L,RANK>;


template <typename L, size_t RANK>
concept liouvillian = hana_sequence<L> && !!hana::all_of(
  decltype(hana::transform(std::declval<L>(), hana::typeid_)){},
  []<class T>(T) { return lindblad<typename T::type,RANK>; });


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


template <size_t RANK>
struct Lindblad
{
  std::string label;
  Jump<RANK> jump;
  Rate<RANK> rate; // could be made std::optional<Rate<RANK>> rate; if it’s null, then the rate is calculated from the jump functional
  Superoperator<RANK> superoperator; // std::optional<Superoperator<RANK>> superoperator; if it’s null, then the superoperator is calculated from the jump functional
};


} // structure
