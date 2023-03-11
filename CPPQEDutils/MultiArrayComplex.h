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
size_t halveExtents(Extents<TWO_TIMES_RANK> extents)
{
  static constexpr size_t RANK=TWO_TIMES_RANK/2;
  return std_ext::ranges::fold(std::views::iota(0ul,RANK),size_t(1), [&] (size_t init, size_t i) {
    size_t a=extents[i];
#ifndef   NDEBUG
    if (size_t b=extents[i+RANK]; a!=b)
      throw std::invalid_argument("MultiArray does not have the shape of a matrix: index no. "+std::to_string(i)+": "+std::to_string(a)+" vs. "+std::to_string(b));
#endif // NDEBUG
    return init*a;
  });
}

template <size_t TWO_TIMES_RANK>
auto matricize(MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = halveExtents(ma.extents);
  return Eigen::Map<CMatrix>{ma.mutableView().dataView.data(),matrixDim,matrixDim};
}

template <size_t TWO_TIMES_RANK>
auto matricize(const MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = halveExtents(ma.extents);
  return Eigen::Map<const CMatrix>{ma.dataView.data(),matrixDim,matrixDim};
}


template <size_t RANK1, size_t RANK2>
auto directProduct(const MultiArray<dcomp,RANK1>& m1, const MultiArray<dcomp,RANK2>& m2,
                   std::function<dcomp(dcomp,dcomp)> func=[](dcomp v1, dcomp v2) {return v1*v2;} )
{
  MultiArray<dcomp,RANK1+RANK2> res{concatenate(m1.extents,m2.extents)};
  for (size_t stride=m1.dataView.size(), i=0; i<stride; ++i) for (size_t j=0; j<m2.dataView.size(); ++j) res.mutableView().dataView[i+stride*j]=func(m1.dataView[i],m2.dataView[j]);
  return res;
}


} // cppqedutils
