// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "StateVector.h"


namespace quantumdata {

  
template <size_t RANK>
using DensityOperatorView=cppqedutils::MultiArrayView<dcomp,2*RANK>;


template <size_t RANK>
using DensityOperatorConstView=cppqedutils::MultiArrayConstView<dcomp,2*RANK>;


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
template<size_t RANK>
struct DensityOperator : ArrayBase<DensityOperator<RANK>>, DimensionsBookkeeper<RANK>
{
  static constexpr size_t N_RANK=RANK;

  using Dimensions=typename DimensionsBookkeeper<RANK>::Dimensions;
  
  typedef ArrayBase<DensityOperator<RANK>> ABase;

  DensityOperator(const DensityOperator&) = delete; DensityOperator& operator=(const DensityOperator&) = delete;
  
  DensityOperator(DensityOperator&&) = default; DensityOperator& operator=(DensityOperator&&) = default;

  DensityOperator(const Dimensions& dimensions, auto&& initializer)
    : ABase{cppqedutils::concatenate(dimensions,dimensions)}, DimensionsBookkeeper<RANK>{this->extents} {initializer(*this);}

  explicit DensityOperator(const Dimensions& dimensions)
    : DensityOperator{dimensions,[](DensityOperator& rho) {rho=0; rho.mutableView().dataView[0]=1.;}} {}
    
  explicit DensityOperator(const StateVector<RANK>& psi) ///< Constructs the class as a dyadic product of `psi` with itself.
    : ABase(psi.dyad()), DimensionsBookkeeper<RANK>{this->extents} {}


  DensityOperator& operator=(const StateVector<RANK>& psi)
  {
    using namespace linalg;
    CMatrix matrix(matrixView());
    CVector vector(psi.vectorView());
    size_t dim(this->getTotalDimension());
    for (size_t i=0; i<dim; i++) for (size_t j=0; j<dim; j++) matrix(i,j)=vector(i)*conj(vector(j));
    return *this;
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

};


template<size_t RANK1, size_t RANK2>
inline
const DensityOperator<RANK1+RANK2>
operator*(const DensityOperator<RANK1>&, const DensityOperator<RANK2>&);


template<size_t RANK>
inline
double frobeniusNorm(const DensityOperator<RANK>& rho) {return rho.frobeniusNorm();}


/// Performs the opposite of quantumdata::deflate
template<size_t RANK>
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
template<size_t RANK>
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


template<typename V, size_t RANK>
auto reduce(const LazyDensityOperator<RANK>& matrix)
{
  return partialTrace<V>(matrix,densityOperatorize<mpl::size<V>::value>);
}


template<size_t... SUBSYSTEM, size_t RANK>
auto reduce(const LazyDensityOperator<RANK>& matrix)
{
  return reduce<tmptools::Vector<SUBSYSTEM...>>(matrix);
}


template<typename V, size_t RANK>
auto reduce(const DensityOperator<RANK>& matrix)
{
  return partialTrace<V>(matrix,[](const auto& v) {return v;});
}


template<size_t... SUBSYSTEM, size_t RANK>
auto reduce(const DensityOperator<RANK>& matrix)
{
  return reduce<tmptools::Vector<SUBSYSTEM...>>(matrix);
}



template<size_t RANK>
inline auto
dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
  return DensityOperator<RANK>(sv1.dyad(sv2),byReference);
}


template <size_t RANK>
constexpr auto MultiArrayRank_v<DensityOperator<RANK>> = 2*RANK;


} // quantumdata


namespace cppqedutils::sliceiterator {

template<size_t RANK>
struct SubscriptMultiArray<quantumdata::DensityOperator<RANK>>
{
  template<typename ... SubscriptPack>
  static auto _(const quantumdata::DensityOperator<RANK>& rho, const SubscriptPack&... subscriptPack)
  {
    return rho.diagonalSliceIndex(subscriptPack...);
  }
};
  
} // cppqedutils::sliceiterator


template<size_t RANK>
double purity(const DensityOperator<RANK>& rho)
{
  return std::accumulate(rho.getArray().begin(),rho.getArray().end(),0.,[] (double init, dcomp element) { return init+cppqedutils::sqrAbs(element); });
}



namespace structure {

using ::quantumdata::DensityOperatorView, ::quantumdata::DensityOperatorConstView;

} // structure

