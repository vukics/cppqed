// -*- C++ -*-
#ifndef   QUANTUMDATA_IMPL_DENSITYOPERATOR_TCC_INCLUDED
#define   QUANTUMDATA_IMPL_DENSITYOPERATOR_TCC_INCLUDED

#include "DensityOperator.h"

#include "impl/BlitzArrayExtensions.tcc"
#include "impl/MultiIndexIterator.tcc"


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
  double res=real(sum(m(i,i)));
  m/=res;
  return res;
}



template<int RANK> 
const dcomp&
DensityOperator<RANK>::operator()(const Idx& i, const Idx& j) const
{
  return operator()()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
  // We have to explicitly indicate the template parameters for concatenateTinies here, otherwise when RANK=1, that is, Idx is int, the compiler cannot find the function. This is despite TinyVector<int,1> has an implicit constructor from int, because implicit type conversions are NEVER considered for template parameter deduction (for further details see EffC++ 3rd edition item 46.)
}



template<int RANK>
void inflate(const TTD_DARRAY(1)& flattened, DensityOperator<RANK>& rho, bool offDiagonals)
{
  using mathutils::sqr;
  typedef typename DensityOperator<RANK>::Dimensions Dimensions;
  
  const size_t dim=rho.getTotalDimension();
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(Dimensions(size_t(0)),rho.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  // Diagonal
  for (Iterator i(etalon); idx<dim; ++i)
    rho(dispatchLDO_index(*i),dispatchLDO_index(*i))=flattened(idx++);
  
  // OffDiagonal
  if (offDiagonals)
    for (Iterator i=etalon.getBegin(); idx<mathutils::sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j, idx+=2) {
        dcomp matrixElement(rho(dispatchLDO_index(*i),dispatchLDO_index(*j))=dcomp(flattened(idx),flattened(idx+1)));
        rho(dispatchLDO_index(*j),dispatchLDO_index(*i))=conj(matrixElement);
      }
  
}


} // quantumdata

#endif // QUANTUMDATA_IMPL_DENSITYOPERATOR_TCC_INCLUDED
