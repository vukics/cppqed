// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDCORE_QUANTUMDATA_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED

#include "NonOrthogonalStateVector.h"

#include "ComplexArrayExtensions.h"

namespace quantumdata {

template<int RANK, typename TRAFO>
double NonOrthogonalStateVector<RANK,TRAFO>::norm() const
{
  return sqrt(real(sum(conj(vectorView())*blitzplusplus::unaryArray(dual_))));
}

#include <iostream>

template<int RANK, typename TRAFO>
void NonOrthogonalStateVector<RANK,TRAFO>::update() const
{
  TrafoTraits::transform(transformation_,operator()(),dual_);
}

} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED

