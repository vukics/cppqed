// -*- C++ -*-
#ifndef QUANTUMDATA_IMPL_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED
#define QUANTUMDATA_IMPL_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED

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

#endif // QUANTUMDATA_IMPL_NONORTHOGONALSTATEVECTOR_TCC_INCLUDED

