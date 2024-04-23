// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "SparseMatrix.h"


void quantumoperator::SparseMatrix::operator()(std::span<const dcomp> psi, std::span<dcomp> dpsidt) const
{
  for (auto&& [j,i,value] : elements) {
#ifndef   NDEBUG
    if (size_t s=std::size(psi); j>=s || i>=s ) std::range_error("SparseMatrix::operator() index values:"+std::to_string(j)+", "+std::to_string(i)+"; extent: "+std::to_string(s));
#endif // NDEBUG
    dpsidt[j]+=value*psi[i];
  }
}

