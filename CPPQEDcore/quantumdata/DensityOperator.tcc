// -*- C++ -*-
#ifndef   QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED
#define   QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED

#include "DensityOperator.h"

#include "BlitzArrayExtensions.tcc"
#include "MultiIndexIterator.tcc"


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
  if (init) *this=0;
}


template<int RANK>
DensityOperator<RANK>::DensityOperator(const DensityOperator& rho) 
  : LDO_Base(rho.getDimensions()),
    ABase(rho.getArray().copy())
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
  return getArray()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
  // We have to explicitly indicate the template parameters for concatenateTinies here, otherwise when RANK=1, that is, Idx is int, the compiler cannot find the function. This is despite TinyVector<int,1> has an implicit constructor from int, because implicit type conversions are NEVER considered for template parameter deduction (for further details see EffC++ 3rd edition item 46.)
}



template<int RANK>
void inflate(const DArray<1>& flattened, DensityOperator<RANK>& rho, bool offDiagonals)
{
  using mathutils::sqr;

  const size_t dim=rho.getTotalDimension();
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(typename DensityOperator<RANK>::Dimensions(size_t(0)),rho.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  typedef typename DensityOperator<RANK>::Idx Idx;
  // Diagonal
  for (Iterator i(etalon); idx<dim; ++i) {
    const Idx ii(dispatchLDO_index(*i));
    rho(ii,ii)=flattened(idx++);
  }
  
  // OffDiagonal
  if (offDiagonals)
    for (Iterator i(etalon); idx<mathutils::sqr(dim); ++i) {
      const Idx ii(dispatchLDO_index(*i));
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j, idx+=2) {
        const Idx jj(dispatchLDO_index(*j));
        dcomp matrixElement(rho(ii,jj)=dcomp(flattened(idx),flattened(idx+1)));
        rho(jj,ii)=conj(matrixElement);
      }
    }
  
}


template<int RANK>
const DensityOperator<RANK>
densityOperatorize(const LazyDensityOperator<RANK>& matrix)
{
  DensityOperator<RANK> res(matrix.getDimension());
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(typename DensityOperator<RANK>::Dimensions(size_t(0)),matrix.getDimensions()-1,cpputils::mii::begin);
  
  typedef typename DensityOperator<RANK>::Idx Idx;
  
  for (Iterator i(etalon); i!=etalon.getEnd(); ++i) {
    const Idx ii(dispatchLDO_index(*i));
    res(ii,ii)=matrix(ii);
    for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
      const Idx jj(dispatchLDO_index(*j));
      dcomp matrixElement(res(ii,jj)=matrix(ii,jj));
      res(jj,ii)=conj(matrixElement);
    }
  }

  return res;
  
}


template<int... SUBSYSTEM, int RANK>
const DensityOperator<mpl::size<tmptools::Vector<SUBSYSTEM...> >::value>
reduce(const LazyDensityOperator<RANK>& matrix)
{
  static const int RES_ARITY=mpl::size<tmptools::Vector<SUBSYSTEM...> >::value;
  return partialTrace<tmptools::Vector<SUBSYSTEM...>,DensityOperator<RES_ARITY> >(matrix,densityOperatorize<RES_ARITY>);
}


} // quantumdata

#endif // QUANTUMDATA_DENSITYOPERATOR_TCC_INCLUDED
