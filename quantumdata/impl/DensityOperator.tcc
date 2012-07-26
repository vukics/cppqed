// -*- C++ -*-
#ifndef   _DENSITY_OPERATOR_IMPL_H
#define   _DENSITY_OPERATOR_IMPL_H

#include "DensityOperator.h"

#include "impl/BlitzTinyExtensions.tcc"

namespace quantumdata {


template<int RANK>
inline
const DensityOperator<RANK>
dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
  return DensityOperator<RANK>(sv1.dyad(sv2),byReference);
}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const StateVector<RANK>& psi) 
  : LDO_Base(psi.getDimensions()),
    ABase(psi.dyad())
{}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const DensityOperatorLow& rho, ByReference) 
  : LDO_Base(blitzplusplus::halfCutTiny(rho.shape())),
    ABase(rho)
{}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const Dimensions& dimensions, bool init)
  : LDO_Base(dimensions),
    ABase(DensityOperatorLow(blitzplusplus::concatenateTinies(dimensions,dimensions)))
{
  if (init) ABase::operator()()=0;
}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const DensityOperator& rho) 
  : LDO_Base(rho.getDimensions()),
    ABase(rho().copy())
{}




template<int RANK>
double
DensityOperator<RANK>::renorm()
{
  using blitz::tensor::i;
  linalg::CMatrix m(matrixView());
  double res=sum(real(m(i,i)));
  m/=res;
  return res;
}



template<int RANK> 
inline
const dcomp 
DensityOperator<RANK>::operator()(const Idx& i, const Idx& j) const
{
  return operator()()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
  // We have to explicitly indicate the template parameters for
  // concatenateTinies here, otherwise when RANK=1, that is, Idx
  // is int, the compiler cannot find the function. This is despite
  // TinyVector<int,1> has an implicit constructor from int, because
  // implicit type conversions are NEVER considered for template
  // parameter deduction (for further details see EffC++ 3rd edition
  // item 46.)
}


} // quantumdata

#endif // _DENSITY_OPERATOR_IMPL_H
