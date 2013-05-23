// -*- C++ -*-
#ifndef   QUANTUMDATA_IMPL_STATEVECTOR_TCC_INCLUDED
#define   QUANTUMDATA_IMPL_STATEVECTOR_TCC_INCLUDED

#include "StateVector.h"

#include "impl/ComplexArrayExtensions.tcc"

namespace quantumdata {


template<int RANK>
StateVector<RANK>::StateVector(const Dimensions& dimensions, bool init)
  : LDO_Base(dimensions),
    ABase(StateVectorLow(dimensions))
{
  if (init) operator()()=0;
}


template<int RANK>
StateVector<RANK>::StateVector(const StateVector& sv) 
  : LDO_Base(sv.getDimensions()), ABase(sv().copy())
{
}


template<int RANK> template<int RANK2>
StateVector<RANK>::StateVector(const StateVector<RANK2>& psi1, const StateVector<RANK-RANK2>& psi2)
  : LDO_Base(blitzplusplus::concatenateTinies(psi1.getDimensions(),psi2.getDimensions())),
    ABase(blitzplusplus::doDirect<blitzplusplus::dodirect::multiplication,RANK2,RANK-RANK2>(psi1(),psi2()))
{
}


template<int RANK>
const typename StateVector<RANK>::DensityOperatorLow 
StateVector<RANK>::dyad(const StateVector& sv) const 
{
  using namespace blitzplusplus;
  return doDirect<dodirect::multiplication,RANK,RANK>(operator()(),StateVectorLow(conj(sv())));
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

  int dim(getTotalDimension());

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


#endif // QUANTUMDATA_IMPL_STATEVECTOR_TCC_INCLUDED
