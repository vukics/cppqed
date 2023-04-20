// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "ArrayBase.h"


namespace quantumdata {

template <size_t RANK> struct DensityOperator; // forward declaration


template <size_t RANK>
using StateVectorView=cppqedutils::MultiArrayView<dcomp,RANK>;


template <size_t RANK>
using StateVectorConstView=cppqedutils::MultiArrayConstView<dcomp,RANK>;


/// State vector of arbitrary arity
/**
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing the initial condition, to be evolved by the quantum dynamics).
 */
template<size_t RANK>
struct StateVector : ArrayBase<StateVector<RANK>>
{
  using ABase = ArrayBase<StateVector<RANK>>;

  using Dimensions=::cppqedutils::Extents<RANK>;

  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;

  StateVector(StateVector&&) = default; StateVector& operator=(StateVector&&) = default;

  StateVector(::cppqedutils::MultiArray<dcomp,RANK>&& ma) : ABase{std::move(ma)} {}

  StateVector(Dimensions dimensions, auto&& initializer) : ABase{dimensions,std::forward<decltype(initializer)>(initializer)} {}

  explicit StateVector(Dimensions dimensions)
    : StateVector{dimensions, [=] () {auto res(::cppqedutils::multiarray::zeroInit<dcomp,RANK>(dimensions)()); res[0]=1.; return res;}} {}

  // StateVector(auto&& ... a) : ArrayBase<StateVector<RANK>>(std::forward<decltype(a)>(a)...) {}

  /// \name Metric
  //@{
    /// Returns the norm \f$\norm\Psi\f$, implemented in terms of ArrayBase::frobeniusNorm.
  friend double norm(const StateVector& sv) {return frobeniusNorm(sv);}
  
  friend double renorm(StateVector& sv) {double res=norm(sv); sv/=res; return res;}
  //@}
  
  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad),
   * which is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF. 
   */
  void addTo(DensityOperator<RANK>& rho, double weight=1.) const;/*
  {
    using namespace linalg;
    CMatrix matrix(rho.matrixView());
    CVector vector(vectorView());
    size_t dim(this->getTotalDimension());
    for (size_t i=0; i<dim; i++) for (size_t j=0; j<dim; j++) matrix(i,j)+=weight*vector(i)*conj(vector(j));
  }*/
  
 
};


/// Forms a dyad with the argument, a rather expensive operation, implemented as a directProduct
template <size_t RANK>
auto dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
#ifndef NDEBUG
  if (sv1.extents != sv2.extents) throw std::runtime_error("Mismatch in StateVector::dyad dimensions");
#endif // NDEBUG
  return cppqedutils::directProduct(sv1, sv2, [] (dcomp v1, dcomp v2) {return v1*conj(v2);} );
}


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& psi1, const StateVector<RANK2>& psi2)
{
  return StateVector<RANK1+RANK2>(cppqedutils::directProduct<RANK1,RANK2>(psi1,psi2));
}


/// Calculates the inner product, relying on StateVector::vectorView
template<size_t RANK>
dcomp braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2);/*
{
  using blitz::tensor::i;
  linalg::CVector temp(psi1.getTotalDimension());
  temp=conj(psi1.vectorView()(i))*psi2.vectorView()(i);
  return sum(temp);
}*/


template <size_t RANK>
constexpr auto multiArrayRank_v<StateVector<RANK>> = RANK;


} // quantumdata


template <size_t RANK>
constexpr auto cppqedutils::passByValue_v<quantumdata::StateVector<RANK>> = false;


namespace structure {

using ::quantumdata::StateVectorView, ::quantumdata::StateVectorConstView;

} // structure
