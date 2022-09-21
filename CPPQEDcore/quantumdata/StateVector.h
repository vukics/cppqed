// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"


namespace quantumdata {


template <size_t RANK>
using StateVectorView=cppqedutils::MultiArrayView<dcomp,RANK>;


template <size_t RANK>
using StateVectorConstView=cppqedutils::MultiArrayConstView<dcomp,RANK>;



/// State vector of arbitrary arity
/**
 * \tparamRANK
 *
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing first the initial condition, to be evolved by the quantum dynamics).
 * 
 */
template<size_t RANK>
struct StateVector : ArrayBase<StateVector<RANK>>, DimensionsBookkeeper<RANK>
{
  static constexpr size_t N_RANK=RANK;

  typedef ArrayBase<StateVector<RANK>> ABase;

  explicit StateVector(const Dimensions& dimensions, std::function<void(StateVector&)> initializer=[](StateVector& psi) {psi=0;})
    : ABase{dimensions}, DimensionsBookkeeper<RANK>{dimensions} {initializer(*this);}

  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;
  
  StateVector(StateVector&&);
  StateVector& operator=(StateVector&&);

  
  StateVector(const StateVector& sv) ///< Copy constructor using by value semantics, that is, deep copy.
    : LDO_Base(sv.getDimensions()), ABase(sv.getArray().copy()) {}

    ///  Move constructor using the original reference semantics of blitz::Array 
    /**
     * \note It doesn’t touch its argument, as there seems to be no efficient way to invalidate that one
     */
  StateVector(StateVector&& sv) : LDO_Base(sv.getDimensions()), ABase(std::move(sv.getArray())) {}

    /// Constructs the class as the direct product of `psi1` and `psi2`, whose arities add up to `RANK`.
    /**
    * The implementation relies on blitzplusplus::concatenateTinies and blitzplusplus::doDirect.
    * \tparam RANK2 the arity of one of the operands
    */
  template<size_t RANK2>
  StateVector(const StateVector<RANK2>& psi1, const StateVector<RANK-RANK2>& psi2)
    : LDO_Base(blitzplusplus::concatenateTinies(psi1.getDimensions(),psi2.getDimensions())),
      ABase(blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication,RANK2,RANK-RANK2>(psi1.getArray(),psi2.getArray())) {}

  StateVector() : LDO_Base(Dimensions{size_t(0)}), ABase() {}
  
  /// Assignment with by-value semantics.
  /** Default assignment doesn't work, because LazyDensityOperator is always purely constant (const DimensionsBookkeeper base). */
  StateVector& operator=(const StateVector&) = default;
 
  StateVector& operator=(StateVector&& sv)
  {
    ABase::operator=(sv.getArray());
    LDO_Base::setDimensions(sv.getDimensions());
    return *this;
  }
 
  /// \name Subscripting
  //@{
  template<typename... SubscriptPack>
  const dcomp& operator()(size_t s0, SubscriptPack... subscriptPack) const ///< Multi-array style subscription. \tparam ...SubscriptPack expected as all integers of number RANK-1 (checked @ compile time)
  {
    static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for StateVector." );
    return getArray()(s0,subscriptPack...);
  }
  
  template<typename... SubscriptPack>
  dcomp& operator()(size_t s0, SubscriptPack... subscriptPack) {return const_cast<dcomp&>(static_cast<const StateVector*>(this)->operator()(s0,subscriptPack...));} ///< ”
  //@}

  /// \name LazyDensityOperator diagonal iteration
  //@{
  template<typename... SubscriptPack>
  auto sliceIndex(const SubscriptPack&... subscriptPack) const
  {
    static_assert( sizeof...(SubscriptPack)==RANK , "Incorrect number of subscripts for StateVector." );
#define SLICE_EXPR getArray()(subscriptPack...)
    return StateVector<cppqedutils::Rank_v<decltype(SLICE_EXPR)>>(SLICE_EXPR,byReference);
#undef  SLICE_EXPR
  }
  
  template<typename... SubscriptPack>
  void transposeSelf(SubscriptPack... subscriptPack)
  {
    static_assert( sizeof...(SubscriptPack)==RANK , "Incorrect number of subscripts for StateVector." );
    getArray().transposeSelf(subscriptPack...);
    this->setDimensions(getArray().shape());
  }
  //@}
  
  /// \name Metric
  //@{
    /// Returns the norm \f$\norm\Psi\f$
    /** Implemented in terms of ArrayBase::frobeniusNorm. */
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
    /// Forms a dyad with the argument
    /** This is a rather expensive operation, implemented in terms of blitzplusplus::doDirect. */
  inline auto dyad(const StateVector& sv) const
  {
    return blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication,RANK,RANK>(getArray(),StateVectorLow(conj(sv.getArray())));
  }

  auto dyad() const {return dyad(*this);} ///< dyad with the object itself
  //@}

  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad).
   * This is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF. 
   */
  void addTo(DensityOperator<RANK>& rho, double weight=1.) const
  {
    using namespace linalg;
    CMatrix matrix(rho.matrixView());
    CVector vector(vectorView());
    size_t dim(this->getTotalDimension());
    for (size_t i=0; i<dim; i++) for (size_t j=0; j<dim; j++) matrix(i,j)+=weight*vector(i)*conj(vector(j));
  }
  
  void reference(const StateVector& other) {getArray().reference(other.getArray()); this->setDimensions(other.getDimensions());}
  
  auto lbound() const {return getArray().lbound();}
  auto ubound() const {return getArray().ubound();}

#ifndef NDEBUG
  void debug() const {std::cerr<<"Debug: "<<getArray()<<std::endl;}
#endif // NDEBUG
  
private:
  const dcomp index(const Idx& i, const Idx& j) const override {return getArray()(i)*conj(getArray()(j));} ///< This function implements the LazyDensityOperator interface in a dyadic-product way.
  
  double trace_v() const override {return norm();} ///< A straightforward implementation of a LazyDensityOperator virtual
  
};


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& t1, const StateVector<RANK2>& t2)
{
  return StateVector<RANK1+RANK2>(t1,t2);
}


/// Calculates the inner product, relying on StateVector::vectorView
template<size_t RANK>
dcomp braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2)
{
  using blitz::tensor::i;
  linalg::CVector temp(psi1.getTotalDimension());
  temp=conj(psi1.vectorView()(i))*psi2.vectorView()(i);
  return sum(temp);
}


template <size_t RANK>
constexpr auto MultiArrayRank_v<StateVector<RANK>> = RANK;


} // quantumdata


namespace structure {

using ::quantumdata::StateVectorView, ::quantumdata::StateVectorConstView;

} // structure
