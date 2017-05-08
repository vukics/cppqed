// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMDATA_STATEVECTOR_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_STATEVECTOR_TCC_INCLUDED

#include "StateVector.h"

#include "DensityOperator.tcc"

#include "ComplexArrayExtensions.tcc"

namespace quantumdata {


template<int RANK>
StateVector<RANK>::StateVector(const Dimensions& dimensions, bool init)
  : LDO_Base(dimensions),
    ABase(StateVectorLow(dimensions))
{
  if (init) *this=0;
}


template<int RANK>
StateVector<RANK>::StateVector(const StateVector& sv) 
  : LDO_Base(sv.getDimensions()), ABase(sv.getArray().copy())
{
}


template<int RANK> template<int RANK2>
StateVector<RANK>::StateVector(const StateVector<RANK2>& psi1, const StateVector<RANK-RANK2>& psi2)
  : LDO_Base(blitzplusplus::concatenateTinies(psi1.getDimensions(),psi2.getDimensions())),
    ABase(blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication,RANK2,RANK-RANK2>(psi1.getArray(),psi2.getArray()))
{
}


template<int RANK> template<typename... SubscriptPack>
const dcomp& StateVector<RANK>::operator()(int s0, SubscriptPack... subscriptPack) const
{
  static_assert( mpl::size<mpl::vector<SubscriptPack...> >::value==RANK-1 , "Incorrect number of subscripts for StateVector." );
  return getArray()(s0,subscriptPack...);
}


template<int RANK>
auto StateVector<RANK>::dyad(const StateVector& sv) const -> const DensityOperatorLow
{
  using namespace blitzplusplus;
  return doDirect<dodirect::multiplication,RANK,RANK>(getArray(),StateVectorLow(conj(sv.getArray())));
}


template<int RANK>
double 
StateVector<RANK>::renorm()
{
  double res=norm();
  operator/=(res); return res;
}


template<int RANK>
void 
StateVector<RANK>::addTo(DensityOperator<RANK>& rho) const
{
  using namespace linalg;

  CMatrix matrix(rho.matrixView());

  const CVector& vector(vectorView());

  int dim(this->getTotalDimension());

  for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) matrix(i,j)+=vector(i)*conj(vector(j));

}


template<int RANK>
const dcomp
braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2)
{
  using blitz::tensor::i;
  linalg::CVector temp(psi1.getTotalDimension());
  temp=conj(psi1.vectorView()(i))*psi2.vectorView()(i);
  return sum(temp);
}


} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_STATEVECTOR_TCC_INCLUDED
