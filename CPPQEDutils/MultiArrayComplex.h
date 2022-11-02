// Copyright András Vukics 2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ComplexExtensions.h"
#include "MultiArray.h"


namespace cppqedutils {


template <size_t RANK1, size_t RANK2>
auto directProduct(const MultiArray<dcomp,RANK1>& m1, const MultiArray<dcomp,RANK2>& m2,
                   std::function<dcomp(dcomp,dcomp)> func=[](dcomp v1, dcomp v2) {return v1*v2;} )
{
  MultiArray<dcomp,RANK1+RANK2> res{concatenate(m1.extents,m2.extents)};
  for (size_t stride=m1.dataView.size(), i=0; i<stride; ++i) for (size_t j=0; j<m2.dataView.size(); ++j) res.mutableView().dataView[i+stride*j]=func(m1.dataView[i],m2.dataView[j]);
  return res;
}


} // cppqedutils
