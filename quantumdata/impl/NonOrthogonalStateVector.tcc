// -*- C++ -*-
#ifndef _NO_STATE_VECTOR_IMPL_H
#define _NO_STATE_VECTOR_IMPL_H

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

#endif // _NO_STATE_VECTOR_IMPL_H

