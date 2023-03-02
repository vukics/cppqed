// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"


namespace quantumdata {

template<size_t RANK> struct DensityOperator; // forward declaration


template <size_t RANK>
using StateVectorView=cppqedutils::MultiArrayView<dcomp,RANK>;


template <typename, size_t RANK> // the first template parameter is there just in order to conform with what MultiArray expects
struct StateVectorConstViewBase : cppqedutils::MultiArrayConstView<dcomp,RANK>
{
  using cppqedutils::MultiArrayConstView<dcomp,RANK>::MultiArrayConstView;
  using cppqedutils::MultiArrayConstView<dcomp,RANK>::operator=;
  
  /// lazy_density_operator-style indexing
  const dcomp operator() (std::convertible_to<size_t> auto ... i) const requires (sizeof...(i)==2*RANK);

};


template <size_t RANK>
using StateVectorConstView = StateVectorConstViewBase<dcomp,RANK>;


/// State vector of arbitrary arity
/**
 * \tparamRANK
 *
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing the initial condition, to be evolved by the quantum dynamics).
 * 
 */
template<size_t RANK>
struct StateVector : ArrayBase<StateVector<RANK>,StateVectorConstViewBase>, DimensionsBookkeeper<RANK>
{
  using Dimensions=typename DimensionsBookkeeper<RANK>::Dimensions;

  using ABase = ArrayBase<StateVector<RANK>,StateVectorConstViewBase>;

  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;
  
  StateVector(StateVector&&) = default; StateVector& operator=(StateVector&&) = default;

  StateVector(const Dimensions& dimensions, auto&& initializer) : ABase{dimensions}, DimensionsBookkeeper<RANK>{dimensions} {initializer(*this);}

  explicit StateVector(const Dimensions& dimensions)
    : StateVector{dimensions,[](StateVector& psi) {psi=0.; psi.mutableView().dataView[0]=1.;}} {}
    
  /// Constructs the class as the direct product of `psi1` and `psi2`, whose arities add up to `RANK`.
  template<size_t RANK2>
  StateVector(const StateVector<RANK2>& psi1, const StateVector<RANK-RANK2>& psi2)
    : ABase{cppqedutils::directProduct(psi1,psi2)}, DimensionsBookkeeper<RANK>{cppqedutils::concatenate(psi1.extents,psi2.extents)} {}
  
  /// \name Metric
  //@{
    /// Returns the norm \f$\norm\Psi\f$, implemented in terms of ArrayBase::frobeniusNorm.
  double   norm() const {return ABase::frobeniusNorm();}
  
  double renorm() ///< ” and also renormalises
  {
    double res=norm();
    this->operator/=(res);
    return res;
  }
  //@}
  
  /// \name Dyad
  //@{
    /// Forms a dyad with the argument, a rather expensive operation, implemented as a directProduct
  inline auto dyad(const StateVector& sv) const
  {
#ifndef NDEBUG
    if (this->getDimensions() != sv.getDimensions()) throw std::runtime_error("Mismatch in StateVector::dyad dimensions");
#endif // NDEBUG
    return cppqedutils::directProduct(*this,sv);
  }

  auto dyad() const {return dyad(*this);} ///< dyad with the object itself
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


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& t1, const StateVector<RANK2>& t2)
{
  return StateVector<RANK1+RANK2>(t1,t2);
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


namespace structure {

using ::quantumdata::StateVectorView, ::quantumdata::StateVectorConstView;

} // structure
