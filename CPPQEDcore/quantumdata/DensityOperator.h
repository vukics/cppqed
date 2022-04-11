// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED

#include "QuantumDataFwd.h"

#include "ArrayBase.h"
#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "BlitzArrayExtensions.h"
#include "MultiIndexIterator.h"


namespace quantumdata {


/// Density operator of arbitrary arity
/**
 * Cf. \ref quantumdatahighlevel "rationale"
 * 
 * \tparamRANK
 * 
 * The DensityOperator interface is similar to StateVector with obvious differences.
 * 
 * \note A DensityOperator `<RANK>` represents a density operator on a Hilbert space of arity `RANK`. This makes that the number of its indices is actually `2*RANK`.
 * 
 */
template<int RANK>
class DensityOperator
  : public LazyDensityOperator<RANK>,
    public ArrayBase<DensityOperator<RANK>>
{
public:
  static const int N_RANK=RANK;
  
  typedef LazyDensityOperator<RANK> LDO_Base;
  
  typedef ArrayBase<DensityOperator<RANK>> ABase;

  typedef typename LDO_Base::Dimensions Dimensions;

  typedef typename LDO_Base::Idx Idx;

  using DensityOperatorLow=typename ABase::ArrayLow;

  typedef linalg::CMatrix CMatrix;

  using ABase::frobeniusNorm; using ABase::getArray; using ABase::operator=;

  /*
  DensityOperator() : Base() {}
  DensityOperatorBase() : LDO_Base(Dimensions()), ABase(DensityOperatorLow()), matrixView_() {}
  */

  DensityOperator(const DensityOperatorLow& rho, ByReference) ///< Referencing constructor implemented in terms of blitzplusplus::halfCutTiny
    : LDO_Base(blitzplusplus::halfCutTiny(rho.shape())), ABase(rho) {}

  explicit DensityOperator(const Dimensions& dimensions, bool init=true)
    : LDO_Base(dimensions),
      ABase(DensityOperatorLow(blitzplusplus::concatenateTinies(dimensions,dimensions))) {if (init) *this=0;}

  explicit DensityOperator(const StateVector<RANK>& psi) ///< Constructs the class as a dyadic product of `psi` with itself.
    : LDO_Base(psi.getDimensions()), ABase(psi.dyad()) {}

  DensityOperator(const DensityOperator& rho) ///< By-value copy constructor (deep copy)
    : LDO_Base(rho.getDimensions()), ABase(rho.getArray().copy()) {}


  DensityOperator(DensityOperator&& rho) ///< Move constructor (shallow copy)
    : LDO_Base(rho.getDimensions()), ABase(std::move(rho.getArray())) {}
    
  DensityOperator() : LDO_Base(Dimensions{size_t(0)}), ABase() {}

  DensityOperator& operator=(const DensityOperator&) = default;

  DensityOperator& operator=(DensityOperator&& rho)
  {
    ABase::operator=(rho.getArray());
    LDO_Base::setDimensions(rho.getDimensions());
    return *this;
  }
  
  DensityOperator& operator=(const StateVector<RANK>& psi)
  {
    using namespace linalg;
    CMatrix matrix(matrixView());
    CVector vector(psi.vectorView());
    int dim(this->getTotalDimension());
    for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) matrix(i,j)=vector(i)*conj(vector(j));
    return *this;
  }
  
private:
  class IndexerProxy
  {
  public:
    IndexerProxy(const DensityOperator* rho, const Idx& firstIndex) : rho_(rho), firstIndex_(firstIndex) {}

    template<typename... SubscriptPack>
    IndexerProxy(const DensityOperator* rho, int s0, SubscriptPack... subscriptPack) : IndexerProxy(rho,Idx(s0,subscriptPack...))
    {static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for DensityOperator." );}

    const dcomp& operator()(const Idx& secondIndex) const {return rho_->indexWithTiny(firstIndex_,secondIndex);}
          dcomp& operator()(const Idx& secondIndex)       {return const_cast<dcomp&>(static_cast<const IndexerProxy&>(*this)(secondIndex));}

    template<typename... SubscriptPack>
    const dcomp& operator()(int s0, SubscriptPack... subscriptPack) const
    {
      static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for DensityOperator." );
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

  /// \name LazyDensityOperator diagonal iteration
  //@{
  template<typename... SubscriptPack>
  auto diagonalSliceIndex(const SubscriptPack&... subscriptPack) const
  {
    static_assert( sizeof...(SubscriptPack)==RANK , "Incorrect number of subscripts for DensityOperator." );
#define SLICE_EXPR getArray()(subscriptPack...,subscriptPack...)
    return DensityOperator<cppqedutils::Rank_v<decltype(SLICE_EXPR)>/2>(SLICE_EXPR,byReference);
#undef  SLICE_EXPR
  }
  
  template<typename... SubscriptPack>
  void transposeSelf(SubscriptPack... subscriptPack)
  {
    static_assert( sizeof...(SubscriptPack)==RANK , "Incorrect number of subscripts for DensityOperator." );
    getArray().transposeSelf(subscriptPack...,(subscriptPack+RANK)...);
    this->setDimensions(blitzplusplus::halfCutTiny(getArray().shape()));
  }
  //@}
  
  /// \name Norm
  //@{
  double norm() const ///< returns the trace norm
  {
    using blitz::tensor::i;
    const linalg::CMatrix m(matrixView());
    return real(sum(m(i,i)));
  }

  double renorm() ///< ” and also renormalises
  {
    double trace=norm();
    *this/=trace;
    return trace;
  }
  //@}

  /// \name Matrix view
  //@{
  auto matrixView() const {return blitzplusplus::binaryArray(getArray());} ///< returns a two-dimensional view of the underlying data, created on the fly via blitzplusplus::binaryArray
  //@}

  void reference(const DensityOperator& other) {getArray().reference(other.getArray()); this->setDimensions(other.getDimensions());}

  auto lbound() const {return blitzplusplus::halfCutTiny(getArray().lbound());}
  auto ubound() const {return blitzplusplus::halfCutTiny(getArray().ubound());}

private:
  const dcomp& indexWithTiny(const Idx& i, const Idx& j) const ///< Used for implementing operator() and the index function below.
  {
    return getArray()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
  }
  
  const dcomp index(const Idx& i, const Idx& j) const override {return indexWithTiny(i,j);} ///< This function implements the LazyDensityOperator interface in a trivial element-access way

  double trace_v() const override {return norm();} ///< A straightforward implementation of a LazyDensityOperator virtual

};


template<int RANK1, int RANK2>
inline
const DensityOperator<RANK1+RANK2>
operator*(const DensityOperator<RANK1>&, const DensityOperator<RANK2>&);


template<int RANK>
inline
double frobeniusNorm(const DensityOperator<RANK>& rho) {return rho.frobeniusNorm();}


/// Performs the opposite of quantumdata::deflate
template<int RANK>
void inflate(const DArray<1>& flattened, DensityOperator<RANK>& rho, bool offDiagonals)
{
  const size_t dim=rho.getTotalDimension();
  
  typedef cppqedutils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(rho.getDimensions()-1,cppqedutils::mii::begin);
  
  size_t idx=0;

  // Diagonal
  for (Iterator i(etalon); idx<dim; ++i)
    rho(*i)(*i)=flattened(idx++);
  
  // OffDiagonal
  if (offDiagonals)
    for (Iterator i(etalon); idx<cppqedutils::sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j, idx+=2) {
        dcomp matrixElement(rho(*i)(*j)=dcomp(flattened(idx),flattened(idx+1)));
        rho(*j)(*i)=conj(matrixElement);
      }

}


/// Creates a DensityOperator as the (deep) copy of the data of a LazyDensityOperator of the same arity
template<int RANK>
DensityOperator<RANK>
densityOperatorize(const LazyDensityOperator<RANK>& matrix)
{
  DensityOperator<RANK> res(matrix.getDimension());
  
  typedef cppqedutils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cppqedutils::mii::begin);
  
  for (Iterator i(etalon); i!=etalon.getEnd(); ++i) {
    res(*i)(*i)=matrix(*i);
    for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
      dcomp matrixElement(res(*i)(*j)=matrix(*i)(*j));
      res(*j)(*i)=conj(matrixElement);
    }
  }

  return res;
  
}


template<typename V, int RANK>
auto reduce(const LazyDensityOperator<RANK>& matrix)
{
  return partialTrace<V>(matrix,densityOperatorize<mpl::size<V>::value>);
}


template<int... SUBSYSTEM, int RANK>
auto reduce(const LazyDensityOperator<RANK>& matrix)
{
  return reduce<tmptools::Vector<SUBSYSTEM...>>(matrix);
}


template<typename V, int RANK>
auto reduce(const DensityOperator<RANK>& matrix)
{
  return partialTrace<V>(matrix,[](const auto& v) {return v;});
}


template<int... SUBSYSTEM, int RANK>
auto reduce(const DensityOperator<RANK>& matrix)
{
  return reduce<tmptools::Vector<SUBSYSTEM...>>(matrix);
}



template<int RANK>
inline auto
dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
  return DensityOperator<RANK>(sv1.dyad(sv2),byReference);
}


template <int RANK>
constexpr auto ArrayRank_v<DensityOperator<RANK>> = 2*RANK;


template<int RANK, typename ... SubscriptPack>
auto subscript(const quantumdata::DensityOperator<RANK>& rho, const SubscriptPack&... subscriptPack) ///< for use in cppqedutils::SliceIterator
{
  return rho.diagonalSliceIndex(subscriptPack...);
}

} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_DENSITYOPERATOR_H_INCLUDED
