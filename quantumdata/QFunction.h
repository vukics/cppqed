// -*- C++ -*-
#ifndef QUANTUMDATA_QFUNCTION_H_INCLUDED
#define QUANTUMDATA_QFUNCTION_H_INCLUDED

#include "BlitzArray.h"
#include "ComplexExtensions.h"
#include "MathExtensions.h"

#include "ParsFwd.h"

namespace quantumdata {


namespace details {

double qFunctionHelper(size_t n, const dcomp& alpha);
  
}


template<typename DensityOperator>
double qFunction(const DensityOperator& rho, double x, double y)
{
  using namespace mathutils;

  const dcomp alpha(x,y);

  dcomp qComplex;
  
  for (size_t m=0; m<rho.getDimension(); ++m) for (size_t n=0; n<rho.getDimension(); ++n)
    qComplex+=details::qFunctionHelper(n,alpha)*details::qFunctionHelper(m,conj(alpha))*rho(m,n);
  
  return exp(sqrAbs(alpha))*real(qComplex);
}
  

}

#endif // QUANTUMDATA_QFUNCTION_H_INCLUDED
