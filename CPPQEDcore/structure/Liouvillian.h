// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "LiouvillianAveragedCommon.h"

#include "LazyDensityOperator.h"
#include "StateVector.h"

#include "VectorFromMatrixSliceIterator.h"

#include <list>
#include <tuple>
#include <variant>


template <int RANK>
using TimeDependentJump = std::function<void(double t, quantumdata::StateVectorLow<RANK>& psi)>;

template <int RANK>
using TimeIndependentJump = std::function<void(quantumdata::StateVectorLow<RANK>& psi)>;


template <int RANK>
using TimeDependentRate = TimeDependentExpectationValue<RANK>;

template <int RANK>
using TimeIndependentRate = TimeIndepedentExpectationValue<RANK>;


template <int RANK>
using TimeDependentSuperoperator = std::function<void(double t, const quantumdata::DensityOperatorLow<RANK>& rho, 
                                                      quantumdata::DensityOperatorLow<RANK>& drhodt)>;

template <int RANK>
using TimeIndependentSuperoperator = std::function<void(const quantumdata::DensityOperatorLow<RANK>& rho,
                                                        quantumdata::DensityOperatorLow<RANK>& drhodt)>;


template <int RANK>
using Rate = std::variant<TimeDependentRate<RANK>,TimeIndependentRate<RANK>>;

template <int RANK>
using Jump = std::variant<TimeDependentJump<RANK>,TimeIndependentJump<RANK>>;

template <int RANK>
using Superoperator = std::variant<TimeDependentSuperoperator<RANK>,TimeIndependentSuperoperator<RANK>>;


template <int RANK>
struct Lindblad
{
  std::string label;
  Jump<RANK> jump;
  std::optional<Rate<RANK>> rate; // if it’s null, then the rate is calculated from the jump functional
  std::optional<Superoperator<RANK>> superoperator; // if it’s null, then the superoperator is calculated from the jump functional
};


template <int RANK>
auto rateFromJump(double t, const quantumdata::StateVector<RANK>& psi, Jump<RANK> jump)
{
  quantumdata::StateVector<RANK> psiTemp(psi);
  if (!jump.index()) std::get<0>(jump)(t,psiTemp.getArray());
  else std::get<1>(jump)(psiTemp.getArray());
  return std::tuple{cppqedutils::sqr(psiTemp.norm()),
                    psiTemp // careful! psiTemp is not normalized
                   };
}


template <int RANK>
void superoperatorFromJump(double t, const quantumdata::DensityOperatorLow<RANK>& rho,
                           quantumdata::DensityOperatorLow<RANK>& drhodt, Jump<RANK> jump)
{
  quantumdata::DensityOperatorLow<RANK> rhotemp(rho.copy());

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
}



template <int RANK>
struct ElementLiouvillian
{
  std::string label;
  std::list<Lindblad<RANK>> lindblads;
};


template <int RANK>
std::ostream& streamKey(std::ostream& os, const ElementLiouvillian<RANK>& el, size_t i) // the implementation of this can be put into a .cc file together with streamKey for ElementAveraged -- can be done with std::variant
{
  using namespace std;
  os<<el.label;
  for (const auto& l : el.lindblads) os<<endl<<setw(2)<<i++<<". "<<l.label;
  return os<<endl;
}

