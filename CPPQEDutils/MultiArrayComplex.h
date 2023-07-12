// Copyright András Vukics 2022–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "MultiArray.h"

#include <Eigen/Dense>


using CMatrix=Eigen::MatrixX<dcomp>;
using CVector=Eigen::VectorX<dcomp>;


namespace cppqedutils {


template <size_t RANK>
auto vectorize(MultiArray<dcomp,RANK>& ma)
{
  return Eigen::Map<CVector>{ma.mutableView().dataView.data(),ma.dataView.size()};
}

template <size_t RANK>
auto vectorize(const MultiArray<dcomp,RANK>& ma)
{
  return Eigen::Map<const CVector>{ma.dataView.data(),ma.dataView.size()};
}


template <size_t TWO_TIMES_RANK> requires ( !(TWO_TIMES_RANK%2) )
Extents<TWO_TIMES_RANK/2> halveExtents(Extents<TWO_TIMES_RANK> extents)
{
  static constexpr size_t RANK=TWO_TIMES_RANK/2;
  Extents<RANK> res;
  std::ranges::for_each(std::views::iota(0ul,RANK), [&] (size_t i) {
    res[i]=extents[i];
#ifndef   NDEBUG
    if (size_t b=extents[i+RANK]; res[i]!=b)
      throw std::invalid_argument("MultiArray does not have the shape of a matrix: index no. "+std::to_string(i)+": "+std::to_string(res[i])+" vs. "+std::to_string(b));
#endif // NDEBUG
  });
  return res;
}

template <size_t TWO_TIMES_RANK>
auto matricize(MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  return Eigen::Map<CMatrix>{ma.mutableView().dataView.data(),matrixDim,matrixDim};
}

template <size_t TWO_TIMES_RANK>
auto matricize(const MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  return Eigen::Map<const CMatrix>{ma.dataView.data(),matrixDim,matrixDim};
}


template <size_t RANK1, size_t RANK2>
auto directProduct(const MultiArray<dcomp,RANK1>& m1, const MultiArray<dcomp,RANK2>& m2,
                   std::function<dcomp(dcomp,dcomp)> func=std::multiplies<dcomp>{} )
{
  MultiArray<dcomp,RANK1+RANK2> res{concatenate(m1.extents,m2.extents)};
  for (size_t stride=m1.dataView.size(), i=0; i<stride; ++i) for (size_t j=0; j<m2.dataView.size(); ++j) res.mutableView().dataView[i+stride*j]=func(m1.dataView[i],m2.dataView[j]);
  return res;
}


template <size_t RANK>
MultiArray<dcomp,RANK>& conj(MultiArray<dcomp,RANK>& ma)
{
  for (dcomp& v : ma.mutableView().dataView) v=conj(v);
  return ma;
}


template <size_t RANK>
double frobeniusNorm(const MultiArray<dcomp,RANK>& ma) { return sqrt( std::ranges::fold_left_first( ma.dataView | std::views::transform(sqrAbs), std::plus{} ).value_or(0.) ); }


template <size_t TWO_TIMES_RANK>
void hermitianConjugateSelf(MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  for (size_t i=0; i<matrixDim; ++i) for (size_t j=i+1; j<matrixDim; ++j) {
    dcomp
      &u=ma.mutableView().dataView[i+matrixDim*j],
      &l=ma.mutableView().dataView[j+matrixDim*i];
    u=conj(u); l=conj(l);
    std::swap(u, l);
  }
}


template <size_t TWO_TIMES_RANK>
void twoTimesRealPartOfSelf(MultiArrayView<dcomp,TWO_TIMES_RANK> mav)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(mav.extents));
  auto _=[&] (size_t i, size_t j) -> dcomp& {return mav.dataView[i+matrixDim*j];};
  for (size_t i=0; i<matrixDim; ++i) {
    _(i,i)=2.*real(_(i,i));
    for (size_t j=i; j<matrixDim; ++j) _(j,i)=conj(_(i,j)=_(i,j)+_(j,i));
  }
};


} // cppqedutils
