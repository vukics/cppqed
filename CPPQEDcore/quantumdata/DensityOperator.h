// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED

#include "DensityOperatorFwd.h"

#include "StateVectorFwd.h"

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "Operators.h"


namespace quantumdata {


template<int RANK>
inline
const DensityOperator<RANK>
dyad(const StateVector<RANK>&, const StateVector<RANK>&);


template<int RANK1, int RANK2>
inline
const DensityOperator<RANK1+RANK2>
operator*(const DensityOperator<RANK1>&, const DensityOperator<RANK2>&);


template<int RANK>
inline
double frobeniusNorm(const DensityOperator<RANK>& rho) {return rho.frobeniusNorm();}


/// Density operator of arbitrary arity
/**
 * Cf. \ref quantumdatahighlevel "rationale"
 * 
 * \tparamRANK
 * 
 * The DensityOperator interface is similar to StateVector with obvious differences.
 * 
 * \note A DensityOperator `<RANK>` represents a density operator on a Hilbert space of arity `RANK`. This makes that
 * the number of its indices is actually `2*RANK`. This is the reason why it inherits from quantumdata::ArrayBase <2*RANK>.
 * 
 * \todo provide a move constructor
 * 
 */
template<int RANK>
class DensityOperator
  : public LazyDensityOperator<RANK>,
    private ArrayBase<2*RANK>,
    private linalg::VectorSpace<DensityOperator<RANK> >
{
public:
  static const int N_RANK=RANK;

  typedef LazyDensityOperator<  RANK> LDO_Base;
  typedef ArrayBase          <2*RANK>    ABase;

  typedef typename LDO_Base::Dimensions Dimensions;

  typedef typename LDO_Base::Idx Idx;

  typedef typename ABase::ArrayLow DensityOperatorLow;

  typedef linalg::CMatrix CMatrix;

  using ABase::frobeniusNorm; using ABase::getArray;

  /*
  DensityOperator() : Base() {}
  DensityOperatorBase() : LDO_Base(Dimensions()), ABase(DensityOperatorLow()), matrixView_() {}
  */

  DensityOperator(const DensityOperatorLow&, ByReference); ///< Referencing constructor implemented in terms of blitzplusplus::halfCutTiny

  explicit DensityOperator(const Dimensions&, bool init=true);

  explicit DensityOperator(const StateVector<RANK>& psi); ///< Constructs the class as a dyadic product of `psi` with itself.

  DensityOperator(const DensityOperator&); ///< By-value copy constructor (deep copy)

  /// Default assignment doesn't work, because LazyDensityOperator is always purely constant (const DimensionsBookkeeper base)
  DensityOperator& operator=(const DensityOperator& rho) {ABase::operator=(rho.getArray()); return *this;}

  template<typename OTHER>
  DensityOperator& operator=(const OTHER& other) {getArray()=other; return *this;}

private:
  class IndexerProxy
  {
  public:
    IndexerProxy(const DensityOperator* rho, const Idx& firstIndex) : rho_(rho), firstIndex_(firstIndex) {}

    template<typename... SubscriptPack>
    IndexerProxy(const DensityOperator* rho, int s0, SubscriptPack... subscriptPack) : IndexerProxy(rho,Idx(s0,subscriptPack...))
    {static_assert( mpl::size<mpl::vector<SubscriptPack...> >::value==RANK-1 , "Incorrect number of subscripts for DensityOperator." );}

    const dcomp& operator()(const Idx& secondIndex) const {return rho_->indexWithTiny(firstIndex_,secondIndex);}
          dcomp& operator()(const Idx& secondIndex)       {return const_cast<dcomp&>(static_cast<const IndexerProxy&>(*this)(secondIndex));}

    template<typename... SubscriptPack>
    const dcomp& operator()(int s0, SubscriptPack... subscriptPack) const
    {
      static_assert( mpl::size<mpl::vector<SubscriptPack...> >::value==RANK-1 , "Incorrect number of subscripts for DensityOperator::IndexerProxy." );
      return operator()(Idx(s0,subscriptPack...));
    }

    template<typename... SubscriptPack>
          dcomp& operator()(int s0, SubscriptPack... subscriptPack) {return const_cast<dcomp&>(static_cast<const IndexerProxy&>(*this)(s0,subscriptPack...));}

  private:
    const DensityOperator*const rho_;
    const Idx firstIndex_;

  };

  friend class IndexerProxy;

public:
  /// \name (Multi-)matrix style indexing
  //@{
  const IndexerProxy operator()(const Idx& firstIndex) const {return IndexerProxy(this,firstIndex);}
        IndexerProxy operator()(const Idx& firstIndex)       {return IndexerProxy(this,firstIndex);}

  template<typename... SubscriptPack>
  const IndexerProxy operator()(int s0, SubscriptPack... subscriptPack) const {return IndexerProxy(this,s0,subscriptPack...);}

  template<typename... SubscriptPack>
        IndexerProxy operator()(int s0, SubscriptPack... subscriptPack)       {return IndexerProxy(this,s0,subscriptPack...);}
  //@}

  /// \name Norm
  //@{
  double   norm() const; ///< returns the *trace* “norm”
  double renorm()      ; ///< ” and also renormalises
  //@}

  /// \name Matrix view
  //@{
  const CMatrix matrixView() const {return blitzplusplus::binaryArray(getArray());} ///< returns a two-dimensional view of the underlying data, created on the fly via blitzplusplus::binaryArray
        CMatrix matrixView()       {return blitzplusplus::binaryArray(getArray());} ///< ”
  //@}
  
  /// \name Naive operations for vector space
  //@{
  DensityOperator& operator+=(const DensityOperator& rho) {ABase::operator+=(rho); return *this;}
  DensityOperator& operator-=(const DensityOperator& rho) {ABase::operator-=(rho); return *this;}

  const DensityOperator operator-() const {DensityOperator res(this->getDimensions(),false); res.getArray()=-this->getArray(); return res;}
  const DensityOperator operator+() const {return *this;}

  template<typename OTHER>
  DensityOperator& operator*=(const OTHER& dc) {ABase::operator*=(dc); return *this;}

  template<typename OTHER>
  DensityOperator& operator/=(const OTHER& dc) {ABase::operator/=(dc); return *this;}
  //@}

private:
  const dcomp& indexWithTiny(const Idx&, const Idx&) const; ///< Used for implementing operator() and the index function below.

  const dcomp index(const Idx& i, const Idx& j) const override {return indexWithTiny(i,j);} ///< This function implements the LazyDensityOperator interface in a trivial element-access way

  double trace_v() const override {return norm();} ///< A straightforward implementation of a LazyDensityOperator virtual

};


/// Performs the opposite of quantumdata::deflate
template<int RANK>
void inflate(const DArray<1>&, DensityOperator<RANK>&, bool offDiagonals);


/// Creates a DensityOperator as the (deep) copy of the data of a LazyDensityOperator of the same arity
template<int RANK>
const DensityOperator<RANK>
densityOperatorize(const LazyDensityOperator<RANK>&);


template<int... SUBSYSTEM, int RANK>
const DensityOperator<mpl::size<tmptools::Vector<SUBSYSTEM...> >::value>
reduce(const LazyDensityOperator<RANK>&);


} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED
